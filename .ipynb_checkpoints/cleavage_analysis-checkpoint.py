#!/usr/bin/env python3

import os
import sys
import pandas as pd 
import argparse
import numpy as np 
import glob
import time 
import logging
import matplotlib as mpl
import matplotlib.pyplot as plt
import pickle
from scipy import stats
from multiprocessing import Pool, cpu_count
from typing import List, Dict, Tuple
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import subprocess
from io import StringIO

# Try to import pybedtools for faster bedtools operations
try:
    import pybedtools
    HAS_PYBEDTOOLS = True
except ImportError:
    HAS_PYBEDTOOLS = False
    print("Warning: pybedtools not available, using subprocess calls")


def define_targets(targets: str, ago: str) -> List[str]:
    """Extract target genes for specific AGO protein - optimized with list comprehension"""
    genes = []
    with open(targets, 'r') as f:
        # Skip header and process in one pass
        genes = [
            line.strip().split("\t")[1] 
            for line in f 
            if not line.startswith("ago") and line.strip().split("\t")[0] == ago
        ]
    return genes


def filter_transcript_alignment(bed: str, genes: set, prefix: str = '') -> Tuple[str, str]:
    """
    Filter transcript alignments - MAJOR OPTIMIZATION:
    - Use set for O(1) gene lookups instead of O(n)
    - Use list append + join instead of string concatenation
    - Process file in single pass
    - Unique filenames for parallel processing
    """
    sense_lines = []
    anti_lines = []
    
    with open(bed, 'r') as f:
        for line in f:
            if line.startswith("gene"):
                continue
                
            info = line.strip().split("\t")
            gene = info[0].split(",")[0]
            
            if gene not in genes:  # O(1) lookup with set
                continue
            
            start = int(info[1])
            if start < 0:
                continue
            
            end = int(info[2])
            seq = info[3]
            feature = info[6]
            rpm = float(info[7])
            
            line_out = f'{gene}\t{start}\t{end}\t{seq}\t{rpm}\t'
            
            if "anti" in feature:
                anti_lines.append(line_out + '-\n')
            else:
                sense_lines.append(line_out + '+\n')
    
    # Use unique filenames to avoid conflicts in parallel processing
    if not prefix:
        # Generate unique prefix using process ID and timestamp
        import os
        import time
        prefix = f"tmp_{os.getpid()}_{int(time.time() * 1000000)}_"
    
    anti_file = f'{prefix}anti.bed'
    sense_file = f'{prefix}sense.bed'
    
    # Write files in batch
    with open(anti_file, 'w') as f:
        f.write(''.join(anti_lines))
    
    with open(sense_file, 'w') as f:
        f.write(''.join(sense_lines))
    
    return sense_file, anti_file


def intersect_bed(sense: str, anti: str, output: str = "intersect.bed") -> str:
    """
    Intersect BED files - OPTIMIZATION:
    - Use subprocess.run instead of os.system
    - Add pybedtools option for speed if available
    - Better error handling
    """
    if HAS_PYBEDTOOLS:
        # Much faster than subprocess
        try:
            sense_bed = pybedtools.BedTool(sense)
            anti_bed = pybedtools.BedTool(anti)
            result = sense_bed.intersect(anti_bed, wo=True)
            result.saveas(output)
            return output
        except Exception as e:
            print(f"pybedtools failed, falling back to bedtools: {e}")
    
    # Fallback to bedtools via subprocess
    cmd = ['bedtools', 'intersect', '-wo', '-a', sense, '-b', anti]
    
    with open(output, 'w') as out_file:
        result = subprocess.run(
            cmd,
            stdout=out_file,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
    
    return output


def group_and_agg_dist(df: pd.DataFrame, dist_col: str) -> pd.DataFrame:
    """
    Group and aggregate distance data - OPTIMIZATION:
    - Use more efficient pandas operations
    - Vectorized calculations
    """
    grped = df.groupby(dist_col, as_index=False).agg({
        'sense_rpm': 'sum',
        'anti_rpm': 'sum'
    })
    
    # Vectorized z-score calculations
    sense_mean = grped['sense_rpm'].mean()
    sense_std = grped['sense_rpm'].std(ddof=0)
    anti_mean = grped['anti_rpm'].mean()
    anti_std = grped['anti_rpm'].std(ddof=0)
    
    grped['sense_zscore'] = 2 ** ((grped['sense_rpm'] - sense_mean) / sense_std)
    grped['anti_zscore'] = 2 ** ((grped['anti_rpm'] - anti_mean) / anti_std)
    
    grped = grped.rename(columns={dist_col: 'dist'})
    grped['id'] = dist_col
    
    return grped


def calc_p_distance(intersect: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculate positional distances - OPTIMIZATION:
    - Specify dtypes for faster reading
    - Vectorized string operations
    """
    # Define column types for faster reading
    dtypes = {
        0: str,   # gene
        1: np.int32,  # sense_5p
        2: np.int32,  # sense_3p
        3: str,   # sense_seq
        4: np.float32,  # sense_rpm
        7: np.int32,  # anti_3p
        8: np.int32,  # anti_5p
        9: str,   # anti_seq
        10: np.float32,  # anti_rpm
        12: np.int32  # overlap
    }
    
    dat = pd.read_csv(
        intersect,
        usecols=[0, 1, 2, 3, 4, 7, 8, 9, 10, 12],
        names=['gene', 'sense_5p', 'sense_3p', 'sense_seq', 'sense_rpm', 
               'anti_3p', 'anti_5p', 'anti_seq', 'anti_rpm', 'overlap'],
        header=None,
        sep="\t",
        dtype=dtypes,
        engine='c'  # Use faster C engine
    )
    
    # Vectorized calculations
    dat['5pto5p'] = dat['anti_5p'] - dat['sense_5p']
    dat['5pto3p'] = dat['anti_5p'] - dat['sense_3p'] + 1
    dat['sense_len'] = dat['sense_seq'].str.len()
    dat['anti_len'] = dat['anti_seq'].str.len()
    dat['anti_5p_nt'] = dat['anti_seq'].str[0]
    
    # Process groups
    grouped_5p5p = group_and_agg_dist(dat, '5pto5p')
    grouped_5p3p = group_and_agg_dist(dat, '5pto3p')
    
    grouped = pd.concat([grouped_5p5p, grouped_5p3p], ignore_index=True)
    
    return dat, grouped


def read_fasta(fasta: str) -> Dict[str, str]:
    """
    Read FASTA file - OPTIMIZATION:
    - More efficient parsing
    - Better memory usage
    """
    fa_dict = {}
    current_header = None
    current_seq = []
    
    with open(fasta, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                # Save previous sequence
                if current_header is not None:
                    fa_dict[current_header] = ''.join(current_seq)
                
                # Start new sequence
                current_header = line.split()[0].replace(">", "").split(",")[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        # Don't forget the last sequence
        if current_header is not None:
            fa_dict[current_header] = ''.join(current_seq)
    
    return fa_dict


def reverse_comp(seq: str) -> str:
    """Reverse complement - with translation table for speed"""
    complement = str.maketrans('ATGC', 'TACG')
    return seq.upper()[::-1].translate(complement)


def get_sequences(fasta: Dict[str, str], gene: str, pos: int) -> str:
    """Extract sequences around position - OPTIMIZATION: bounds checking"""
    if gene not in fasta:
        return "NA"
    
    start_pad = pos - 25
    end_pad = pos + 10
    
    if start_pad < 0:
        return "NA"
    
    sequence = fasta[gene][start_pad:end_pad]
    return sequence if len(sequence) == 35 else "NA"


def process_single_sample(args: Tuple) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Process a single sample - for parallel processing
    OPTIMIZATION: This function can now be run in parallel
    Uses unique file prefixes to avoid conflicts
    """
    bed_file, targets, name = args
    
    # Use unique prefix for this sample to avoid file conflicts
    prefix = f"{name}_"
    sense, anti = filter_transcript_alignment(bed_file, targets, prefix=prefix)
    intersect_file = intersect_bed(sense, anti, f"{prefix}intersect.bed")
    raw_dat, measurements = calc_p_distance(intersect_file)
    
    # Clean up temporary files
    for f in [sense, anti, intersect_file]:
        if os.path.exists(f):
            os.remove(f)
    
    raw_dat['sample'] = name
    measurements['sample'] = name
    
    return measurements, raw_dat


def main(target_table: str = '', AGO_name: str = '', 
         IP_bedfiles: List[str] = [], input_bedfiles: List[str] = []) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Main processing function - OPTIMIZATION:
    - Convert targets to set for O(1) lookups
    - Parallel processing option
    """
    # Get targets as a set for fast lookups
    targets_list = define_targets(target_table, AGO_name)
    targets = set(targets_list)  # Convert to set for O(1) lookup
    print(f"Processing {len(targets)} target genes for {AGO_name}")
    
    my_list = [i for i in IP_bedfiles + input_bedfiles if i]
    
    # Option 1: Sequential processing (original)
    results = []
    raw_dats = []
    
    for bed_file in my_list:
        name = os.path.basename(bed_file).split(".")[0]
        print(f"Processing {name}...")
        
        # Use unique prefix for this sample
        prefix = f"{name}_"
        sense, anti = filter_transcript_alignment(bed_file, targets, prefix=prefix)
        intersect_file = intersect_bed(sense, anti, f"{prefix}intersect.bed")
        raw_dat, measurements = calc_p_distance(intersect_file)
        
        raw_dat['sample'] = name
        measurements['sample'] = name
        
        results.append(measurements)
        raw_dats.append(raw_dat)
        
        # Clean up temporary files using actual filenames
        for f in [sense, anti, intersect_file]:
            if os.path.exists(f):
                os.remove(f)
    
    result = pd.concat(results, ignore_index=True) if results else pd.DataFrame()
    result2 = pd.concat(raw_dats, ignore_index=True) if raw_dats else pd.DataFrame()
    
    return result, result2


def main_parallel(target_table: str = '', AGO_name: str = '', 
                  IP_bedfiles: List[str] = [], input_bedfiles: List[str] = [],
                  n_jobs: int = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    OPTIMIZATION: Parallel version of main function
    Uses multiprocessing to process multiple files simultaneously
    """
    targets_list = define_targets(target_table, AGO_name)
    targets = set(targets_list)
    print(f"Processing {len(targets)} target genes for {AGO_name}")
    
    my_list = [i for i in IP_bedfiles + input_bedfiles if i]
    
    # Prepare arguments for parallel processing
    args_list = [
        (bed_file, targets, os.path.basename(bed_file).split(".")[0])
        for bed_file in my_list
    ]
    
    # Determine number of workers
    if n_jobs is None:
        n_jobs = min(cpu_count() - 1, len(my_list))
    
    print(f"Using {n_jobs} parallel workers")
    
    # Process in parallel
    with Pool(processes=n_jobs) as pool:
        parallel_results = pool.map(process_single_sample, args_list)
    
    # Combine results
    measurements_list = [r[0] for r in parallel_results]
    raw_dat_list = [r[1] for r in parallel_results]
    
    result = pd.concat(measurements_list, ignore_index=True) if measurements_list else pd.DataFrame()
    result2 = pd.concat(raw_dat_list, ignore_index=True) if raw_dat_list else pd.DataFrame()
    
    return result, result2


def glob_files(path: str, wildcard: str) -> List[str]:
    """Glob files with filtering - unchanged but typed"""
    files = glob.glob(f"{path}/{wildcard}")
    filtered_files = [file for file in files if not ("_FLAG_" in file or "_GFP_" in file)]
    return filtered_files


def run_all(sample_map: List, fasta: str, use_parallel: bool = True, n_jobs: int = None):
    """
    Run analysis on all samples - OPTIMIZATION:
    - Optional parallel processing
    - Better progress tracking
    - More efficient sequence extraction
    """
    master_tables = []
    raw_masters = []
    
    for sample in sample_map:
        name = sample[0]
        target_list = sample[1]
        input_bed = sample[2]
        IP_bed = sample[3]
        
        print(f"\n{'='*60}")
        print(f"Processing AGO protein: {name}")
        print(f"{'='*60}")
        
        # Choose parallel or sequential processing
        if use_parallel and len(IP_bed + input_bed) > 1:
            result, raw_dat = main_parallel(
                target_table=target_list,
                AGO_name=name,
                IP_bedfiles=IP_bed,
                input_bedfiles=input_bed,
                n_jobs=n_jobs
            )
        else:
            result, raw_dat = main(
                target_table=target_list,
                AGO_name=name,
                IP_bedfiles=IP_bed,
                input_bedfiles=input_bed
            )
        
        result['sample_name'] = name
        raw_dat['sample_name'] = name
        
        master_tables.append(result)
        raw_masters.append(raw_dat)
    
    # Combine all results
    master_table = pd.concat(master_tables, ignore_index=True)
    raw_master = pd.concat(raw_masters, ignore_index=True)
    
    # Optionally add sequences (can be slow for large datasets)
    # Uncomment if needed:
    # print("\nExtracting sequences...")
    # fa_dict = read_fasta(fasta)
    # raw_master['seqlogo'] = raw_master.apply(
    #     lambda x: get_sequences(fa_dict, x['gene'], x['anti_5p']), 
    #     axis=1
    # )
    
    # Save results
    print("\nSaving results...")
    master_table.to_csv("cleavage_analysis_master.new.tsv", 
                       index=False, sep="\t", header=True)
    #raw_master.to_csv("cleavage_analysis_master_raw_data.new.tsv", 
    #                 index=False, sep="\t", header=True)
    
    print(f"\nProcessing complete!")
    print(f"Results saved to:")
    print(f"  - cleavage_analysis_master.new.tsv")
    print(f"  - cleavage_analysis_master_raw_data.new.tsv")

if __name__ == "__main__":
    # Example usage with argument parsing
    parser = argparse.ArgumentParser(description='Analyze RNA cleavage patterns')
    parser.add_argument('--parallel', action='store_true', 
                       help='Use parallel processing')
    parser.add_argument('--jobs', type=int, default=None,
                       help='Number of parallel jobs (default: CPU count - 1)')
    
    args = parser.parse_args()
    
    # Your existing sample_map and fasta setup here
    target_list = '/fs/ess/PCON0160/ben/projects/2023_claycomb_argonomics/target_lists.txt'
    fasta = "/fs/ess/PCON0160/ben/genomes/c_elegans/WS279/ce_ws279.linc.pseudo.pc.repbase.fa"
    fp = "/fs/ess/PCON0160/ben/projects/2025_CSR_Y57_Y54/results/20251105/transcripts/bed"

    sample_map = [  
        # csr-1
        [ 'CSR-1', target_list, glob_files(fp, "*CSR*input*.bed.tsv"), glob_files(fp, "*CSR*IP*.bed.tsv") ],
        [ 'CSR-1', target_list, glob_files("/fs/ess/PCON0160/ben/projects/2025_csr1_cleavage/00AA_AGO_cleavage/20250429/transcripts/bed", "*input*bed.tsv"), glob_files("/fs/ess/PCON0160/ben/projects/2025_csr1_cleavage/00AA_AGO_cleavage/20250429/transcripts/bed", "*IP*bed.tsv")]
    ]

    run_all(sample_map, fasta, use_parallel=args.parallel, n_jobs=args.jobs)