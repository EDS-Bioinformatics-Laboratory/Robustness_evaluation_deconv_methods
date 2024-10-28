################################################################################
# SPOTlight implementation
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

# number of cells per cell type to be consider
# 0 if all cells are used
n_cells <- 0


# reading the ST datasets
spatial.count.matrix <- readRDS(paste0(Spatial_data, "count.matrix.", z, ".rds"))
spatial.count.matrix <- exprs(Biobase::ExpressionSet(assayData = as.matrix(spatial.count.matrix)))
spots.metadata <- readRDS(paste0(Spatial_data, "metadata.", z, ".rds"))

# number of spots in the simulated data
n.spots <- dim(spots.metadata)[1]

start_time <- Sys.time()
SPOTlight.method <- function(sc.ref.data, n_cells) {
  # converting Seurat object to SingleCellExperiment
  sc.ref.data.sce <- as.SingleCellExperiment(sc.ref.data)
  colLabels(sc.ref.data.sce) <- colData(sc.ref.data.sce)$blue.main
  
  # create dataframe for marker genes
  mgs <- scoreMarkers(sc.ref.data.sce)
  mgs_fil <- lapply(names(mgs), function(i) {
    x <- mgs[[i]]
    # Filter and keep relevant marker genes, those with AUC > 0.5
    # (changed from 0.8 as in originally suggested by developer)
    x <- x[x$mean.AUC > 0.5,]
    # Sort the genes from highest to lowest weight
    x <- x[order(x$mean.AUC, decreasing = TRUE),]
    # Add gene and cluster id to the dataframe
    x$gene <- rownames(x)
    x$cluster <- i
    data.frame(x)
  })
  mgs_df <- do.call(rbind, mgs_fil)
  # rm(sc.ref.data.sce, mgs_df, mgs_fil)
  
  if (n_cells > 0) {
    # downsampling number of cells per celltype
    idx <- split(seq(ncol(sc.ref.data)), sc.ref.data$blue.main)
    # n_cells <- 20
    cs.keep <- lapply(idx, function(i) {
      n <- length(i)
      if (n < n_cells) {
        n_cells <- n
      }
      sample(i, n_cells)
    })
    sc.ref.data.spotlight <- sc.ref.data[, unlist(cs.keep)]
  }
  else {
    sc.ref.data.spotlight <- sc.ref.data
  }
  
  # running SPOTlight deconvolution method
  message("Running SPOTlight deconvolution algorithm")
  results.SPOTlight <- SPOTlight::SPOTlight(
    x = sc.ref.data.spotlight,
    # single-cell reference data counts
    y = spatial.count.matrix,
    # spatial data count matrix
    groups = sc.ref.data.spotlight$blue.main,
    mgs = mgs_df,
    weight_id = "mean.AUC"
  )
  return(results.SPOTlight)
}
print("Starting SPOTlight deconvolution")
results.SPOTlight <- SPOTlight.method(sc.ref.data, n_cells)
end_time <- Sys.time()
print(paste0(
  "/n Time taken for running SPOTlight algorithm: ",
  round(difftime(end_time, start_time, units = "mins"), 2),
  " min"
))
rm(start_time, end_time)


results.SPOTlight <- as.data.frame(results.SPOTlight$mat)

write.csv(results.SPOTlight,
          file = paste0(Method_Res, "SPOTlight/results.SPOTlight.", z, "-", y, ".csv"))