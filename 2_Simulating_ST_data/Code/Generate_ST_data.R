#### Header start ##############################################################
# We create simulated ST datasets using the single cell RNA-sequencing reference
# data with 13 cell types (cell clusters) and at max 300 cells per cell type
# 
# This script generates spatial count matrix data, spatial spot metadata
# information, and spatial co-ordinate metadata
# 
# Count matrix data: genes * spot matrix
# Spatial spot metadata: spots * cell types matrix
# Spatial co-ordinate metadata: spots * X-Y coordinate columns
#
#
# command line arguments
# variable parameters for generation of ST datasets
# 1. min number of cell types per spot
# 2. max number of cell types per spot
# 3. min number of cells per spot
# 4. max number of cells per spot
# 5. index of simulated ST dataset (either 1/2/3)
#
# Implemented by Utkarsh Mahamune
# Bioinformatics Laboratory | Amsterdam UMC (Location AMC)
#### Header end ################################################################


#### Initialize the environment ####
source("Init_env.R")

start_time <- Sys.time()

args <- commandArgs(trailingOnly = TRUE)
# 1. min number of cell types per spot
# 2. max number of cell types per spot
# 3. min number of cells per spot
# 4. max number of cells per spot
# 5. index of simulated ST dataset (either 1/2/3)


# fixed parameters for generation of ST datasets
# 1. path to single-cell reference data
# 2. path to save spatial transcritomics data + spot metadata + co-odinate metadata
# 3. number of spots to simulate for the spatial transcriptomics data
fix_args <- c("../../1_Generate_sc_ref_data/Results/",
              "../Results/Spatial_Data/",
              1600)


#### scRNA-seq dataset ####
sc.ref.data <- readRDS(paste0(fix_args[1], "sc.ref.data.rds"))
Idents(sc.ref.data) <- "blue.main"

# cell types in the single-cell reference data without NA in the list
celltypes <- sort(unique(sc.ref.data$blue.main))

#### Simulated spatial transcriptomics data ####
# number of spots in the simulated data
n.spots <- fix_args[3]


# for reproducibility
set.seed(167)

# collect all cells of a cell type to generate count matrix and
# randomly select 5 cells out of it for calculating average read counts
get_counts <- function(celltypes, sc.ref.data) {
  ## Reference count matrix for all cell types in the single-cell reference data
  count.celltypes <- list()
  
  for (i in 1:length(celltypes)) {
    # gathering all cells of a celltype
    df <- sc.ref.data@meta.data %>% data.frame()
    cell.names <- df[df$blue.main %in% c(celltypes[i]),] %>% rownames()
    
    # randomly choosing 5 cells from the above compiled list
    random.cell.names <- sample(cell.names, 5)
    
    # count matrix for randomly selected 5 cells of the cell type in consideration
    count_ <-
      as.matrix(sc.ref.data@assays$RNA@counts[, random.cell.names])
    
    # average counts based on 5 randomly selected cells
    counts_ <- ceiling(rowSums(as.matrix(rowSums(count_))) / 5)
    count.celltypes[[i]] <- counts_
  }
  
  count.celltypes <- ((count.celltypes)) %>%
    base::Reduce(function(m1, m2) {
      cbind(unlist(m1), unlist(m2))
    }, .) %>%
    data.frame()
  
  colnames(count.celltypes) <- celltypes
  
  return(count.celltypes)
}


# get spot distribution and metadata
generate_spot_distribution <- lapply(1:n.spots, function(i) {
  # Select number of cell types randomly from the single-cell reference data
  # between min and max cell types per spot
  celltype.pool <- sample(celltypes, sample(x = c(args[1]:args[2]), size = 1))
  
  # probability of each cell type
  prob_ <- stats::runif(length(celltype.pool), 0, 1)
  prob_ <- prob_ / sum(prob_) # normalized sum
  
  # number of cells in the spot between min and max cells per spot
  nr_cells <- sample(args[3]:args[4], 1)
  
  if (nr_cells < length(celltype.pool)) {
    nr_cells <- length(celltype.pool)
  }
  
  # number of cells of a cell type
  nr_cell_per_spot <- c(0)
  while (any(nr_cell_per_spot <= 0)) {
    # prob has been normalized
    nr_cell_per_spot <- stats::rmultinom(n = 1, size = nr_cells, prob = prob_)
  }
  
  # probabilities to get the number of cells per celltype by multinomial 
  # distribution is used as spots metadata
  spot.metadata.prob <- t(as.matrix(prob_)) %>%
    as.data.frame(row.names = paste0("Spot_", i), check.names = T)
  colnames(spot.metadata.prob) <- celltype.pool
  
  # proportions for each celltype present in a spot after multinomial distribution
  # is used as spots metadata
  spot.metadata <- t(as.matrix(nr_cell_per_spot/sum(nr_cell_per_spot))) %>%
    as.data.frame(row.names = paste0("Spot_", i), check.names = T)
  
  colnames(spot.metadata) <- celltype.pool
  
  
  # getting average counts for a cell type using 5 random cells of the same
  # celltype, we aim to have different cells for each spot and thus calculate
  # counts for every spot instead of making one dataframe for all spot
  count.celltypes <- get_counts(celltype.pool, sc.ref.data)
  
  spot.count <- 0
  for (k in 1:length(celltype.pool)) {
    spot.count <- spot.count +
      (rowSums(as.matrix(count.celltypes[celltype.pool[k]]) * nr_cell_per_spot[k]))
  }
  
  # if observed UMIs in spatial transcriptomics data are greater than 30k,
  # downsample counts to random number between 25000 - 30000
  if (sum(spot.count) > 30000) {
    spot.count <- scuttle::downsampleMatrix(
      Matrix::Matrix(spot.count),
      prop = sample(25000:30000, 1) / sum(spot.count),
      bycol = T
    )
  } else {
    spot.count <- Matrix::Matrix(spot.count)
  }
  
  rownames(spot.count) <- names(rowSums(as.matrix(rowSums(count.celltypes))))
  colnames(spot.count) <- paste0("Spot_", i)
  
  return(list(spot.count, spot.metadata, spot.metadata.prob))
})

## getting the count matrix for all spots
spatial.count.matrix <-
  map(generate_spot_distribution, 1) %>% base::Reduce(function(m1, m2) {
    cbind(unlist(m1), unlist(m2))
  }, .)

## getting spot metadata for all spots with proportions from multinomial dist
spots.metadata <- map(generate_spot_distribution, 2) %>%
  dplyr::bind_rows() %>%
  data.frame(check.names = F)

# check for NA values
spots.metadata[is.na(spots.metadata)] <- 0

# adding all cell types columns for the spot metadata
if (sum(celltypes %in% colnames(spots.metadata)) == length(celltypes)) {
  spots.metadata <- spots.metadata[, celltypes]
} else {
  missing_labels <-
    celltypes[which(!celltypes %in% colnames(spots.metadata))]
  # if a cell type is absent set the value to 0
  spots.metadata[missing_labels] <- 0
  spots.metadata <- spots.metadata[, celltypes]
}

## getting spot metadata for all spots with probabilities used for multinomial dist
spots.metadata.prob <- map(generate_spot_distribution, 3) %>%
  dplyr::bind_rows() %>%
  data.frame(check.names = F)

# check for NA values
spots.metadata.prob[is.na(spots.metadata.prob)] <- 0

# adding all cell types columns for the spot metadata
if (sum(celltypes %in% colnames(spots.metadata.prob)) == length(celltypes)) {
  spots.metadata.prob <- spots.metadata.prob[, celltypes]
} else {
  missing_labels <-
    celltypes[which(!celltypes %in% colnames(spots.metadata.prob))]
  # if a cell type is absent set the value to 0
  spots.metadata.prob[missing_labels] <- 0
  spots.metadata.prob <- spots.metadata.prob[, celltypes]
}


## get roughly equal number of rows and columns for the ST data co-ordinates
n.spots <- as.integer(n.spots)
row_col <- function(n.spots) {
  c = as.integer(n.spots ** 0.5 + 0.5)
  while((n.spots %% c) != 0) {
    c = c-1
  }
  return(list(c, n.spots/c))
}
r_c <- unlist(row_col(n.spots))

# spatial transcriptomics co-ordinates
x.coords <- (sample(x = r_c[1]*2, size = n.spots, replace = T))
y.coords <- (sample(x = r_c[2]*2, size = n.spots, replace = T))

spatial.coord <- list(x = x.coords, y = y.coords) %>% data.frame()
rownames(spatial.coord) <- colnames(spatial.count.matrix)

end_time <- Sys.time()
print(message("Time taken for creating simulated ST dataset: ",
              round(difftime(end_time, start_time, units = "mins"), 2), " min"))


# saving the simulated ST datasets components
saveRDS(spatial.count.matrix, file = paste0(fix_args[2], "count.matrix.", args[5], ".rds"))
saveRDS(spots.metadata, file = paste0(fix_args[2], "metadata.", args[5], ".rds"))
saveRDS(spots.metadata.prob, file = paste0(fix_args[2], "metadata.prob.", args[5], ".rds"))
saveRDS(spatial.coord, file = paste0(fix_args[2], "coord.", args[5], ".rds"))

# generating .h5ad version of spatial (count matrix) data
spatial.count.matrix.h5 <- CreateSeuratObject(counts = (spatial.count.matrix),
                                              assay = "Spatial")

SeuratDisk::SaveH5Seurat(spatial.count.matrix.h5, overwrite = T,
             filename = paste0(fix_args[2], "spatial_data.", args[5], ".h5Seurat"))

SeuratDisk::Convert(paste0(fix_args[2], "spatial_data.", args[5], ".h5Seurat"),
        dest = "h5ad",
        overwrite = T)



