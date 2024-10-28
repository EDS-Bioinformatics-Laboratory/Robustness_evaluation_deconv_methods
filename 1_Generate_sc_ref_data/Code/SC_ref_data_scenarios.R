#### Header start ##############################################################
# In the analysis, we require different versions of scRNA-seq reference datasets
# to evaluate the robustness of the deconvolution method
# 
# These versions differ by the cell types present in the data.
# The script generates all such versions using the complete scRNA-seq reference
# data generated after the processing via SC_ref_data.R script.
#
# Implemented by Utkarsh Mahamune
# Bioinformatics Laboratory | Amsterdam UMC (Location AMC)
#### Header end ################################################################

#### Initialize the environment ####
source("Init_env.R")

# generating new sc ref data after removal of cell type/s
new_sc_data <- function(sc.ref.data_, ident_, invert_) {
  sc.ref.data.subset <- subset(sc.ref.data_, idents = ident_, invert = invert_)

  sc.ref.data.subset <- DietSeurat(sc.ref.data.subset)
  sc.ref.data.subset <- ScaleData(sc.ref.data.subset, verbose = F,
                                  features = rownames(sc.ref.data.subset))
  sc.ref.data.subset <- RunPCA(sc.ref.data.subset,
                               features = VariableFeatures(object = sc.ref.data.subset),
                               verbose = F) %>%
    FindNeighbors(dims = 1:30, verbose = F) %>%
    FindClusters(resolution = 0.5, verbose = F) %>%
    RunUMAP(dims = 1:30, seed.use = 15, verbose = F)

  return(sc.ref.data.subset)
}

start_time <- Sys.time()

## original scRNA-seq dataset
org.sc.ref.data <- readRDS(paste0(Results, "sc.ref.data.rds"))


print("Entire single-cell reference data")
Idents(org.sc.ref.data) <- "blue.main"

# creating duplicate of the original scRNA-seq reference data with different
# name so as to fit in the downstream analysis
saveRDS(org.sc.ref.data, file = paste0(Results, "rm0/sc.ref.data.1.rds"))
print(org.sc.ref.data)
sceasy::convertFormat(org.sc.ref.data,
                      from = "seurat",
                      to = "anndata",
                      main_layer = "counts",
                      outFile = paste0(Results, "rm0/sc.ref.data.1.h5ad"),
                      drop_single_values = F)
print("converted")


# cell types in the scRNA-seq reference dataset without NA in the list
celltypes <- sort(unique(org.sc.ref.data$blue.main))

## scRNA-seq reference datasets with one cell type missing
# path to data directory
sc_data <- paste0(Results, "rm1/")
Idents(org.sc.ref.data) <- "blue.main"

for (i in 1:length(celltypes)) {
  sc.ref.data <- new_sc_data(org.sc.ref.data, ident_ = celltypes[i], invert_ = TRUE)
  Idents(sc.ref.data) <- "blue.main"
  saveRDS(sc.ref.data, file = paste0(sc_data, "sc.ref.data.", i, ".rds"))
  sceasy::convertFormat(
    sc.ref.data,
    from = "seurat",
    to = "anndata",
    main_layer = "counts",
    outFile = paste0(sc_data, "sc.ref.data.", i, ".h5ad"),
    drop_single_values = F
  )
}


## scRNA-seq reference datasets with two cell types missing
# path to data directory
sc_data <- paste0(Results, "rm2/")

idents.to.remove <- list(
  c("CD4_T_cells", "CD8_T_cells"), c("Myocytes", "Skeletal_muscle"),
  c("HSC", "NK_cells"), c("Macrophages", "Monocytes"), c("Adipocytes", "Fibroblasts")
)

for (i in 1:length(idents.to.remove)) {
  sc.ref.data_ <- new_sc_data(org.sc.ref.data, ident_ = unlist(idents.to.remove[i]), invert_ = TRUE)
  Idents(sc.ref.data_) <- "blue.main"
  saveRDS(sc.ref.data_, file = paste0(sc_data, "sc.ref.data.", i, ".rds"))
  sceasy::convertFormat(
    sc.ref.data_,
    from = "seurat",
    to = "anndata",
    main_layer = "counts",
    outFile = paste0(sc_data, "sc.ref.data.", i, ".h5ad"),
    drop_single_values = F
  )
}

## scRNA-seq reference datasets with three cell types missing
# path to data directory
sc_data <- paste0(Results, "rm3/")

idents.to.remove <- list(
  c("NK_cells", "CD4_T_cells", "CD8_T_cells"), c("Myocytes", "Skeletal_muscle", "Fibroblasts"),
  c("B_cells", "CD4_T_cells", "CD8_T_cells"), c("Macrophages", "Monocytes", "NK_cells"),
  c("Adipocytes", "Fibroblasts", "Endothelial_cells")
)

for (i in 1:length(idents.to.remove)) {
  sc.ref.data_ <- new_sc_data(org.sc.ref.data, ident_ = unlist(idents.to.remove[i]), invert_ = TRUE)
  Idents(sc.ref.data_) <- "blue.main"
  saveRDS(sc.ref.data_, file = paste0(sc_data, "sc.ref.data.", i, ".rds"))
  sceasy::convertFormat(sc.ref.data_,
                        from = "seurat",
                        to = "anndata",
                        main_layer = "counts",
                        outFile = paste0(sc_data, "sc.ref.data.", i, ".h5ad"),
                        drop_single_values = F)
}

## scRNA-seq reference datasets with five cell types missing
# path to data directory
sc_data <- paste0(Results, "rm5/")

idents.to.remove <- list(
  c("NK_cells", "CD4_T_cells", "CD8_T_cells", "HSC", "B_cells")
)

for (i in 1:length(idents.to.remove)) {
  sc.ref.data_ <- new_sc_data(org.sc.ref.data, ident_ = unlist(idents.to.remove[i]), invert_ = TRUE)
  Idents(sc.ref.data_) <- "blue.main"
  saveRDS(sc.ref.data_, file = paste0(sc_data, "sc.ref.data.", i, ".rds"))
  
  sceasy::convertFormat(sc.ref.data_,
                        from = "seurat",
                        to = "anndata",
                        main_layer = "counts",
                        outFile = paste0(sc_data, "sc.ref.data.", i, ".h5ad"),
                        drop_single_values = F)
}

## scRNA-seq reference datasets with ten cell types missing
# path to data directory
sc_data <- paste0(Results, "rm10/")

idents.to.remove <- list(
  c("NK_cells", "CD4_T_cells", "CD8_T_cells"), c("Myocytes", "Skeletal_muscle", "Fibroblasts"),
  c("B_cells", "CD4_T_cells", "CD8_T_cells"), c("Macrophages", "Monocytes", "NK_cells"),
  c("Adipocytes", "Fibroblasts", "Endothelial_cells")
)

for (i in 1:length(idents.to.remove)) {
  sc.ref.data_ <- new_sc_data(org.sc.ref.data, ident_ = unlist(idents.to.remove[i]), invert_ = FALSE)
  Idents(sc.ref.data_) <- "blue.main"
  saveRDS(sc.ref.data_, file = paste0(sc_data, "sc.ref.data.", i, ".rds"))
  
  sceasy::convertFormat(sc.ref.data_,
                        from = "seurat",
                        to = "anndata",
                        main_layer = "counts",
                        outFile = paste0(sc_data, "sc.ref.data.", i, ".h5ad"),
                        drop_single_values = F)
}

## scRNA-seq reference datasets with eleven cell types missing
# path to data directory
sc_data <- paste0(Results, "rm11/")

idents.to.remove <- list(
  c("CD4_T_cells", "CD8_T_cells"), c("Myocytes", "Skeletal_muscle"),
  c("HSC", "NK_cells"), c("Macrophages", "Monocytes"), c("Adipocytes", "Fibroblasts")
)

for (i in 1:length(idents.to.remove)) {
  sc.ref.data_ <- new_sc_data(org.sc.ref.data, ident_ = unlist(idents.to.remove[i]), invert_ = FALSE)
  Idents(sc.ref.data_) <- "blue.main"
  saveRDS(sc.ref.data_, file = paste0(sc_data, "sc.ref.data.", i, ".rds"))
  sceasy::convertFormat(sc.ref.data_,
                        from = "seurat",
                        to = "anndata",
                        main_layer = "counts",
                        outFile = paste0(sc_data, "sc.ref.data.", i, ".h5ad"),
                        drop_single_values = F)
}

end_time <- Sys.time()
print(paste0(
  "/n Time taken for generating single cell reference dataset version: ",
  round(difftime(end_time, start_time, units = "mins"), 2),
  " min"
))