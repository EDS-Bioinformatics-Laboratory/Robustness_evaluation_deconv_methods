################################################################################
# MuSiC implementation
# 
# Requirement:
# BiocManager::install("TOAST")
# 
# References:
# https://github.com/xuranw/MuSiC
# https://xuranw.github.io/MuSiC/articles/MuSiC.html#reference
# https://www.biorxiv.org/content/10.1101/2022.05.08.491077v1.full
# 
# Implemented by Utkarsh Mahamune
# Bioinformatics Laboratory | Amsterdam UMC (Location AMC)
################################################################################

#### Initialize the environment ####
source("Init_env.R")

# index of simulated ST dataset
# index of single-cell reference dataset
# path to st data
# path to sc reference data
# path to saving results

args <- commandArgs(trailingOnly = TRUE)
z <- args[1]
y <- args[2]
Spatial_data <- args[3]
sc_data <- args[4]
Method_Res <- args[5]

#### scRNA-seq dataset ####
sc.ref.data <- readRDS(paste0(sc_data, "sc.ref.data.", y, ".rds"))
Idents(sc.ref.data) <- "blue.main"

# cell types in the single-cell reference data without NA in the list
celltypes <- unique(sc.ref.data$blue.main)


spatial.count.matrix <- readRDS(paste0(Spatial_data, "count.matrix.", z, ".rds"))
spots.metadata <- readRDS(paste0(Spatial_data, "metadata.", z, ".rds"))

# number of spots in the simulated data
n.spots <- dim(spots.metadata)[1]

start_time <- Sys.time()

# adding cells names as samples under cells column in single cell reference data
sc.ref.data@meta.data$cells <- rownames(sc.ref.data@meta.data)

exp_sc <- as.SingleCellExperiment(sc.ref.data)

exp_st <- exprs(Biobase::ExpressionSet(assayData = as.matrix(spatial.count.matrix)))

results.MuSiC <- MuSiC::music_prop(bulk.mtx = exp_st,
                                   sc.sce = exp_sc,
                                   markers = NULL,
                                   clusters = "blue.main",
                                   samples = "cells",
                                   select.ct = NULL,
                                   cell_size = NULL,
                                   ct.cov = FALSE,
                                   verbose = TRUE)

end_time <- Sys.time()

print(paste0(
  "/n Time taken for running MuSiC algorithm: ",
  round(difftime(end_time, start_time, units = "mins"), 2),
  " min"
))
rm(start_time, end_time)

results.MuSiC <- as.matrix(results.MuSiC$Est.prop.weighted)
results.MuSiC <- results.MuSiC[, order(colnames(results.MuSiC))]

write.csv(results.MuSiC,
          file = paste0(Method_Res, "MuSiC/results.MuSiC.", z, "-", y, ".csv"))