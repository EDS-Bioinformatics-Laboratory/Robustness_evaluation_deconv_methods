  ################################################################################
# Seurat implementation
#
# Referred from
# https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium
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

spatial.count.matrix <-
  readRDS(paste0(Spatial_data, "count.matrix.", z, ".rds"))
spots.metadata <-
  readRDS(paste0(Spatial_data, "metadata.", z, ".rds"))

# number of spots in the simulated data
n.spots <- dim(spots.metadata)[1]

start_time <- Sys.time()

spatial.seurat <-
  CreateSeuratObject(spatial.count.matrix, assay = "Spatial")


spatial.seurat <-
  NormalizeData(spatial.seurat, assay = "Spatial", verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000, verbose = F) %>%
  ScaleData(verbose = F) %>%
  RunPCA(verbose = FALSE)

sc.ref.data.sct <- 
  NormalizeData(sc.ref.data, verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000, verbose = F) %>%
  ScaleData(verbose = F) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30, verbose = FALSE)

anchors <-
  Seurat::FindTransferAnchors(
    reference = sc.ref.data.sct,
    query = spatial.seurat,
    verbose = FALSE
  )

predictions.assay <-
  Seurat::TransferData(
    anchorset = anchors,
    refdata = sc.ref.data.sct$blue.main,
    prediction.assay = TRUE,
    weight.reduction = spatial.seurat[["pca"]],
    dims = 1:30,
    verbose = FALSE
  )
spatial.seurat[["predictions"]] <- predictions.assay

end_time <- Sys.time()
print(paste0(
  "Time taken for running Seurat algorithm: ",
  round(difftime(end_time, start_time, units = "mins"), 2),
  " min"
))

results.Seurat <- as.data.frame(t(spatial.seurat@assays$predictions@data))
# dropping the last column named 'max' from the results which gives rowSum
results.Seurat <- results.Seurat[,-dim(results.Seurat)[2]]
results.Seurat <- results.Seurat[, order(colnames(results.Seurat))]

write.csv(results.Seurat,
          file = paste0(Method_Res, "Seurat/results.Seurat.", z, "-", y, ".csv"))