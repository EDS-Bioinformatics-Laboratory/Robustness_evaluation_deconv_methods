################################################################################
# SCDC implementation
# 
# Rerefences:
# https://github.com/meichendong/SCDC
# 
# The implementation uses SCDC_prop_one() function to find the cell types in a 
# sample (here we have spot since we are dealing with spatial transcriptomics
# data)
# We use weight.basis as FALSE to avoid "Not enough valid cell types" error
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

args <- c(1, 1, "../../2_Simulating_ST_data/Results/Spatial_Data/",
          "../../1_Generate_sc_ref_data/Results/",
          "../Results/")

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

# reading spatial transcriptomics data
spatial.count.matrix <- readRDS(paste0(Spatial_data, "count.matrix.", z, ".rds"))
spots.metadata <- readRDS(paste0(Spatial_data, "metadata.", z, ".rds"))

# number of spots in the simulated data
n.spots <- dim(spots.metadata)[1]

start_time <- Sys.time()

print("Starting SCDC  deconvolution")
# converting scRNa-seq reference data to ExpressionSet type
exp_sc <- Biobase::ExpressionSet(assayData = as.matrix(sc.ref.data@assays$RNA@counts),
                                 phenoData = Biobase::AnnotatedDataFrame(sc.ref.data@meta.data))
# converting simulated spatial transcriptomics data to ExpressionSet type
exp_st <- Biobase::ExpressionSet(assayData = as.matrix(spatial.count.matrix))







bulk.eset <- exp_st[rowSums(exprs(exp_st)) > 0, , drop = FALSE]

sc.basis <- SCDC::SCDC_basis_ONE(x = exp_sc, ct.sub = unique(sc.ref.data$blue.main), 
               ct.varname = "blue.main", sample = "cells", ct.cell.size = NULL)

if (weight.basis) {
  basis <- sc.basis$basis.mvw
} else {
  basis <- sc.basis$basis
}
commongenes <- intersect(rownames(basis), rownames(bulk.eset))
if (length(commongenes) < 0.2 * min(dim(sc.eset)[1], dim(bulk.eset)[1])) {
  stop("Too few common genes!")
}
message(paste("Used", length(commongenes), "common genes..."))
basis.mvw <- basis[commongenes, ct.sub]
xbulk <- getCPM0(exprs(bulk.eset)[commongenes, ])
ALS.S <- sc.basis$sum.mat[ct.sub]
N.bulk <- ncol(bulk.eset)
valid.ct <- (colSums(is.na(basis.mvw)) == 0) & (!is.na(ALS.S))

if (sum(valid.ct) <= 1) {
  stop("Not enough valid cell type!")
}












# The default `weight.basis = TRUE` has changed to FALSE since it gives
# "Not enough valid cell type!" error with weight basis
results.SCDC <- SCDC::SCDC_prop_ONE(bulk.eset = exp_st,
                                    sc.eset = exp_sc,
                                    ct.varname = "blue.main",
                                    sample = "cells",
                                    ct.sub = unique(sc.ref.data$blue.main),
                                    weight.basis = F)


end_time <- Sys.time()
print(paste0(
  "/n Time taken for running SCDC algorithm: ",
  round(difftime(end_time, start_time, units = "mins"), 2),
  " min"
))

results.SCDC <- as.matrix(results.SCDC$prop.est.mvw)
results.SCDC <- results.SCDC[, order(colnames(results.SCDC))]

write.csv(results.SCDC,
          file = paste0(Method_Res, "SCDC/results.SCDC.", z, "-", y, ".csv"))
