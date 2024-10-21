################################################################################
# CARD (Conditional Autoregressive-based Deconvolution) implementation
#
# References
# Publication: https://www.nature.com/articles/s41587-022-01273-7#code-availability
# Tutorial: https://yingma0107.github.io/CARD/documentation/04_CARD_Example.html
#
# Implemented by Utkarsh Mahamune
# Bioinformatics Laboratory | Amsterdam UMC (Location AMC)
################################################################################

#### Initialize the environment ####
source("Init_env.R")

library("wrMisc")

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
spatial.meta.data <- readRDS(paste0(Spatial_data, "coord.", z, ".rds"))

# number of spots in the simulated data
n.spots <- dim(spots.metadata)[1]

start_time <- Sys.time()

# single cell reference data meta data
sc.ref.data.metadata <-
  sc.ref.data@meta.data %>% dplyr::select(blue.main) %>% data.frame()
sc.ref.data.metadata$sampleInfo <- "sample1" # suggested by developer of CARD
head(sc.ref.data.metadata)

print("before CARD obejct")
# creating card object
CARD_obj = CARD::createCARDObject(
  sc_count = sc.ref.data@assays$RNA@counts,
  sc_meta = sc.ref.data.metadata,
  spatial_count = spatial.count.matrix,
  spatial_location = spatial.meta.data,
  ct.varname = "blue.main",
  ct.select = unique(sc.ref.data$blue.main),
  sample.varname = "sampleInfo",
  minCountGene = 100,
  minCountSpot = 5
)

message("Deconvolution starting")
# deconvolution using CARD
results.CARD <- CARD::CARD_deconvolution(CARD_object = CARD_obj)

end_time <- Sys.time()
print(paste0(
  "Time taken for running CARD algorithm: ",
  round(difftime(end_time, start_time, units = "mins"), 2), " min"))


results.CARD <- as.matrix(results.CARD@Proportion_CARD)
results.CARD <- results.CARD[, order(colnames(results.CARD))]

print(getwd())

write.csv(results.CARD,
          file = paste0(Method_Res, "CARD/results.CARD.", z, "-", y, ".csv"))