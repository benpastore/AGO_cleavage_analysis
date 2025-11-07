##############################################
# Install and load required R packages
##############################################

# List of required packages
packages <- c(
  "ggplot2",
  "RColorBrewer",
  "ggrepel",
  "dplyr",
  "scales",
  "ggpubr",
  "tidyr",
  "readr",
  "rstatix",
  "gghalves",
  "ggbeeswarm"
)

# Function to install missing packages
install_if_missing <- function(pkg_list) {
  missing_pkgs <- pkg_list[!(pkg_list %in% installed.packages()[,"Package"])]
  if (length(missing_pkgs)) {
    message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
    install.packages(missing_pkgs, dependencies = TRUE, repos = "https://cloud.r-project.org/")
  } else {
    message("All packages already installed.")
  }
}

# Install and load packages
install_if_missing(packages)

# Load all packages
invisible(lapply(packages, library, character.only = TRUE))

message("✅ All required packages are installed and loaded successfully!")
