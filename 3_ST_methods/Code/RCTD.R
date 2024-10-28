################################################################################
# RCTD implementation
#
# First part also includes creating count matrix for single cell reference data
# which is obsolete in the current version of the scripts and thus commented
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

# library("doSNOW")
# library("doParallel")

renv_lib_paths <- readLines("../../renv_library_paths.txt")
Sys.setenv("R_LIBS_USER"=renv_lib_paths[1])
print(Sys.getenv("R_LIBS_USER"))

args <- commandArgs(trailingOnly = TRUE)
z <- args[1]
y <- args[2]
Spatial_data <- args[3]
sc_data <- args[4]
Method_Res <- args[5]


#### scRNA-seq dataset ####
sc.ref.data <- readRDS(paste0(sc_data, "sc.ref.data.", y, ".rds"))
Idents(sc.ref.data) <- "blue.main"

# # RCTD works with raw counts, thus we use RNA assay from the integrated dataset
# # by subsetting the dataset on the genes and the cells list 
# sc.ref.data.raw <- readRDS("../../1_Generate_sc_ref_data/Results/sc.combined.nor.rds")
# sc.ref.data.raw <- subset(sc.ref.data.raw,
#                           features = rownames(sc.ref.data@assays$RNA@counts),
#                           cells = colnames(sc.ref.data@assays$RNA@counts))
# 
# sc.ref.data <- sc.ref.data.raw

# cell types in the single-cell reference data without NA in the list
celltypes <- unique(sc.ref.data$blue.main)

spatial.count.matrix <- readRDS(paste0(Spatial_data, "count.matrix.", z, ".rds"))
spots.metadata <- readRDS(paste0(Spatial_data, "metadata.", z, ".rds"))
spatial.meta.data <- readRDS(paste0(Spatial_data, "coord.", z, ".rds"))

# number of spots in the simulated data
n.spots <- dim(spots.metadata)[1]

start_time <- Sys.time()

# single cell metadata
barcodes <- Cells(sc.ref.data)
sc.data.meta.data <-
  sc.ref.data@meta.data %>%
  dplyr::select(blue.main, nCount_RNA) %>%
  data.frame(row.names = NULL)
sc.data.meta.data$barcode <- barcodes
sc.data.meta.data <- sc.data.meta.data[, c(3, 1, 2)]
colnames(sc.data.meta.data) = c("barcode", "cluster", "nUMI")

# reference constructor
print("Creating reference constructor")
cell_types.RCTD <- sc.data.meta.data$cluster
names(cell_types.RCTD) <- sc.data.meta.data$barcode
cell_types.RCTD <- as.factor(cell_types.RCTD)
nUMI.RCTD <- sc.data.meta.data$nUMI
names(nUMI.RCTD) <- sc.data.meta.data$barcode

reference.RCTD <- spacexr::Reference(counts = sc.ref.data@assays$RNA@counts,
                                     cell_types = cell_types.RCTD,
                                     nUMI = nUMI.RCTD,
                                     n_max_cells = 10000)

# adding barcodes column to the spatial coordinate metadata
spatial.meta.data <- cbind(rownames(spots.metadata), spatial.meta.data)
colnames(spatial.meta.data) <- c("barcodes", "xcoord", "ycoord")

# spatialRNA constructor
print("Creating spatial RNA constructor")
rownames(spatial.meta.data) <- spatial.meta.data$barcodes
spatial.meta.data$barcodes <- NULL

puck.RCTD <- spacexr::SpatialRNA(coords = spatial.meta.data,
                                 counts = spatial.count.matrix,
                                 nUMI = colSums(spatial.count.matrix))

# creating RCTD object and executing RCTD
print("Running RCTD(spacexr) deconvolution algorithm")
object.RCTD <-
  create.RCTD(puck.RCTD, reference.RCTD)
results.RCTD <- run.RCTD(object.RCTD, doublet_mode = "full")
end_time <- Sys.time()
print(paste0(
  "/n Time taken for running RCTD algorithm: ",
  round(difftime(end_time, start_time, units = "mins"), 2),
  " min"
))

results.RCTD <- as.data.frame(normalize_weights(results.RCTD@results$weights))

write.csv(results.RCTD,
          file = paste0(Method_Res, "RCTD/results.RCTD.", z, "-", y, ".csv"))