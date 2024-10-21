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

# Setting up current working directory with path of this file
if (rstudioapi::isAvailable()) {
  message("running in RStudio: current directory set as working directory")
  currentPath = rstudioapi::getActiveDocumentContext()$path
  setwd(dirname(currentPath))
} else {
  message("running in Terminal")
}


# Paths to different processing folders for better accessibility
Data <- "../../Data/Processed/"
Results <- "../Results/"
Settings <- "../Settings/"


# creating folders for scRNA-seq reference datasets with different scenarios for
# cell type removal to be used during the project
if (dir.exists(Results) != T) {
  dir.create(Results, mode = "0777")
}
if (dir.exists(paste0(Results, "rm1")) != T) {
  dir.create(paste0(Results, "rm1"), mode = "0777")
  dir.create(paste0(Results, "rm2"), mode = "0777")
  dir.create(paste0(Results, "rm3"), mode = "0777")
  dir.create(paste0(Results, "rm5"), mode = "0777")
  dir.create(paste0(Results, "rm10"), mode = "0777")
  dir.create(paste0(Results, "rm11"), mode = "0777")
}


## TO DO ## 
# install packages from different sources if not already present



# Loading and citing the packages
# List of packages required throughout the directory
packagesList <-
  c(
    "ggplot2", # to beautify base R plottings - 3.3.6
    # "fmsb", # to plot ranking/spider/radar plots - 0.7.5
    "tidyverse", # - 1.3.1
    "abind", # - 1.4.5
    # "envalysis", # ggplot theme for publication worthy plots - 0.5.4
    # "dplyr", # no need to install separately - 1.0.9
    # "tibble", # no need to install separately - 3.1.8
    # "tidyr", # no need to install separately - 1.2.0
    "Matrix", # no need to install separately - 1.4-1
    "philentropy", # to calculate JSD b/w probability distributions - 0.6.0
    "Metrics", # to calcualte RMSE b/w probability distributions - 0.1.4
    "viridis", # colors in the plots - 0.6.2
    "SingleCellExperiment", # - 1.18.0
    "SingleR", # annotations package - 1.10.0
    "gtools", # - 3.9.2.1
    # "alluvial", # alluvial/sankey plots - 0.1-2
    "ggalluvial", # 0.12.5
    # "ggrepel", # 0.9.3
    "scran", # functions for interpretation of scRNA-seq - 1.24.0
    # "zellkonverter", # for H5AD file generation - 1.6.1
    "SPOTlight", # deconvolution method - 1.0.0
    "spacexr", # RCTD # deconvolution method, previously known as RCTD - 2.0.0
    "SeuratObject",
    "Seurat", # SeuratObject needs to be installed as a prerequisite  # deconvolution method - 4.1.0
    "SeuratDisk", # dependent on Seurat package - 0.0.0.9019
    "CARD", # deconvolution method - 1.0
    # "wrMisc", # necessary for CARD algorithm - 1.9.1
    "RColorBrewer", # for color pallette - 1.1-3
    # "reticulate", # running python modules in R - 1.25
    "cowplot", # for plotting multiple plots in one figure - 1.1.1
    "gridExtra", # for arranging plots in grid - 2.3
    "grid", # for arranging plots in grid - 4.2.0
    "caret", # for min-max normalization function
    "SCDC", # deconvolution method for bulk RNA-seq data - 0.0.0.9000
    "MuSiC", # deconvolution method for bulk RNA-seq data - 1.0.0
    # "xbioc", # extra utility functions to perform common tasks in analysis - 0.1.19
    "funkyheatmap", # funky plots for summarizing all the results
    "filesstrings", # to move files
    "clustree", # cluster maps for different resolutions
    "sceasy" # converting seurat object to anndata object - 0.0.6
  )


# write citations for the packages loaded in the session
fileConn <- file(paste0(Settings, "citation_renv_R4.1.2.txt"), "wt")

for (pkg in packagesList) {
  # library(pkg, character.only = TRUE)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  # suppressWarnings(writeLines(c(capture.output(toBibtex(citation(pkg)))), fileConn))
  message("Package ", pkg, " with version ", packageVersion(pkg))
}
close(fileConn)

# write session information and citations for each run
writeLines(
  capture.output(sessionInfo()),
  paste0(Settings, "sessionInfo_renv_R4.1.2.txt")
)

# Removing unnecessary variables
rm(pkg, fileConn, packagesList)