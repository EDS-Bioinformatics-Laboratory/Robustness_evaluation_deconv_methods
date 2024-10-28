################################################################################
# Initialization file common to all the scripts in the directory
# This script should be included in all the other scripts
#
# loads R libraries, setup paths to different directories
#
# Implemented by Utkarsh Mahamune
# Bioinformatics Laboratory | Amsterdam UMC (Location AMC)
################################################################################

# Clearing work space
rm(list = ls())

# adding renv library paths to .libPaths()
renv_lib_paths <- readLines("../../renv_library_paths.txt")
.libPaths(c(.libPaths(), renv_lib_paths))

# if (rstudioapi::isAvailable()) {
#   message("running in RStudio: current directory set as working directory")
#   currentPath = rstudioapi::getActiveDocumentContext()$path
#   setwd(dirname(currentPath))
# }

# Paths to different processing folders for better accessibility
Data <- "../../Data/Processed/"
Results <- "../Results/"
Settings <- "../Settings/"
logs <- "logs"

if (dir.exists(Results) != T) {
  dir.create(Results, mode = "0777")
}

if (dir.exists(logs) != T) {
  dir.create(logs, mode = "0777")
}

# Loading and citing the packages
# List of packages required throughout the directory
packagesList <-
  c(
    "ggplot2", # to beautify base R plottings - 3.3.6
    "tidyverse", # - 1.3.1
    "abind", # - 1.4.5
    "Matrix", # no need to install separately - 1.4-1
    "philentropy", # to calculate JSD b/w probability distributions - 0.6.0
    "Metrics", # to calcualte RMSE b/w probability distributions - 0.1.4
    "viridis", # colors in the plots - 0.6.2
    "SingleCellExperiment", # - 1.18.0
    "SingleR", # annotations package - 1.10.0
    "gtools", # - 3.9.2.1
    "ggalluvial", # 0.12.5
    "scran", # functions for interpretation of scRNA-seq - 1.24.0
    "SPOTlight", # deconvolution method - 1.0.0
    "spacexr", # RCTD # deconvolution method, previously known as RCTD - 2.0.0
    "SeuratObject",
    "Seurat", # SeuratObject needs to be installed as a prerequisite  # deconvolution method - 4.1.0
    "SeuratDisk", # dependent on Seurat package - 0.0.0.9019
    "CARD", # deconvolution method - 1.0
    "RColorBrewer", # for color pallette - 1.1-3
    "cowplot", # for plotting multiple plots in one figure - 1.1.1
    "gridExtra", # for arranging plots in grid - 2.3
    "grid", # for arranging plots in grid - 4.2.0
    "caret", # for min-max normalization function
    "SCDC", # deconvolution method for bulk RNA-seq data - 0.0.0.9000
    "MuSiC", # deconvolution method for bulk RNA-seq data - 1.0.0
    "funkyheatmap", # funky plots for summarizing all the results
    "filesstrings", # to move files
    "clustree", # cluster maps for different resolutions
    "sceasy" # converting seurat object to anndata object - 0.0.6
  )


for (pkg in packagesList) {
  suppressWarnings(suppressPackageStartupMessages(library(pkg, character.only = TRUE)))
}

# write session information for each run
sink(paste0(Settings, "sessionInfo.txt"))
sessionInfo()
sink()


# Removing unnecessary variables
rm(pkg, packagesList)