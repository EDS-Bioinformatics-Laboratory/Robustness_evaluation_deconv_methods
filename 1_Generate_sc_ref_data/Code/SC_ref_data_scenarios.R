#### Header start ##############################################################
# In the analysis, we require different versions of scRNA-seq reference data to
# evaluate the robustness of the deconvolution methods.
# These versions differ by the cell types present in the data.
# The script generates all such versions using the complete scRNA-seq reference
# data generated after the processing via SC_ref_data_donorwise.R script.
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
    FindNeighbors(dims = 1:50, verbose = F) %>%
    FindClusters(resolution = 0.5, verbose = F) %>%
    RunUMAP(dims = 1:50, seed.use = 42, verbose = F)

  return(sc.ref.data.subset)
}


## original scRNA-seq dataset
org.sc.ref.data <- readRDS(paste0(Results, "sc.ref.data.rds"))

# org.sc.ref.data <- readRDS(paste0(Results, "sc.combined.nor.rds"))


# sc.ref.data <- org.sc.ref.data

print("Entire single-cell reference data")
Idents(org.sc.ref.data) <- "blue.main"

# creating duplicate of the original scRNA-seq reference data with different
# name so as to fit in the downstream analysis
if (file.exists(paste0(Results, "sc.ref.data.1.rds")) != T) {
  saveRDS(org.sc.ref.data, file = paste0(Results, "sc.ref.data.1.rds"))
}

# creating .h5ad format version of the scRNA-seq reference dataset
if (file.exists(paste0(Results, "sc.ref.data.1.h5ad")) != T) {
  sceasy::convertFormat(org.sc.ref.data,
                        from = "seurat",
                        to = "anndata",
                        main_layer = "counts",
                        outFile = paste0(Results, "sc.ref.data.1.h5ad"),
                        drop_single_values = F)
}


# cell types in the scRNA-seq reference dataset without NA in the list
celltypes <- sort(unique(org.sc.ref.data$blue.main))

## scRNA-seq reference datasets with one cell type missing
# path to data directory
sc_data <- paste0(Results, "rm1/")
Idents(org.sc.ref.data) <- "blue.main"

print("Removing 1 celltype")
for (i in 1:length(celltypes)) {
  sc.ref.data <- new_sc_data(org.sc.ref.data, ident_ = celltypes[i], invert_ = TRUE)
  Idents(sc.ref.data) <- "blue.main"
  print(sort(unique(sc.ref.data$blue.main)))
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

print("2")
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

print("3")
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


print("5")
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

print("10")
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

print("11")
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