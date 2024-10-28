#### Header start ##############################################################
# The script splits the tabula sapiens data donor-wise and platform-wise and
# consider them as an independent scRNA-seq datasets.
# 
# Along the line, the script will save R objects in `../Results` and they are
# treated as checkpoints to resume analysis
# Various plots and figures are saved in the `../Results` directory as well
# This script will also generate UMAP representation at various stages of data
# processing and bar plot visualization for understanding the composition of data
# The plots can be found in '../Results' directory as well.
# 
# 
# References, quality control:
# https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html
# 
# Implemented by Utkarsh Mahamune
# Bioinformatics Laboratory | Amsterdam UMC (Location AMC)
#### Header end ################################################################



#### Initialize the environment ####
source("Init_env.R")

#### Cluster annotation function ####
# annotate the clusters using SingleR R package
annotate_clusters <- function(sc.object) {
  # Load the HumanCellAtlas database
  HumanCellAtlas <- celldex::HumanPrimaryCellAtlasData()
  Blueprint <- celldex::BlueprintEncodeData()
  
  # For the individual databases converting the Seurat object to a
  # 'SingleCellExperiment' format
  seurat.sce <- Seurat::as.SingleCellExperiment(sc.object, assay = "RNA")
  
  dbList <- list(human = HumanCellAtlas, blue = Blueprint)
  
  for (grain in c("main")) {
    for (i in 1:length(dbList)) {
      db <- dbList[i]
      # And extract the features we are going to use for cell identification
      common <- intersect(rownames(seurat.sce), rownames(db[[1]]))
      seurat.db <- seurat.sce[common,]
      
      # Test and reference sets should always be log normalized.
      # The included reference sets are already normalized.
      seurat.db <- scater::logNormCounts(seurat.db)
      
      pred.db <- SingleR::SingleR(test = seurat.db,
                                  ref = db[[1]][common, ],
                                  labels = db[[1]][common, ]$label.main,
                                  assay.type.ref = "logcounts")
      
      # Annotate the original seurat object with the labels
      sc.object@meta.data[, paste0(names(db), '.', grain)] <- pred.db$pruned.labels
    }
  }
  return(sc.object)
}

#### Color pallette function ####
get_pallete <- function(n) {
  # https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
  # Maximum distinct colors you can get is 28
  library(RColorBrewer)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual' &
                                    brewer.pal.info$colorblind == T,]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors,
                             rownames(qual_col_pals)))
  return(col_vector[1:n])
}


#### Processing data function ####
# the function performs several steps of processing functionality in a pipeline
process_data <- function(unprocessed.data) {
  processed.data <- NormalizeData(unprocessed.data, verbose = F) %>%
    FindVariableFeatures(selection.method = "vst",
                         nfeatures = 3000,
                         verbose = F) %>%
    ScaleData(verbose = F) %>%
    RunPCA(npcs = 30, verbose = F) %>%
    RunUMAP(reduction = "pca", dims = 1:30, verbose = F)
  
  return(processed.data)
}


#### Tabula sapiens data ####

# Reading tabula sapiens data
if (file.exists(paste0(Results, "TS_Lymph_Node.h5seurat")) != T) {
  # converting .h5ad file to .h5seurat format to use with seurat package in R
  message("Converting data file from '.h5ad' to '.h5seurat' version")
  SeuratDisk::Convert(paste0(Data, "TS_Lymph_Node.h5ad"), dest = "h5seurat")
  
  # moving the resulted file to Results directory of the particular computation
  filesstrings::file.move(paste0(Data, "TS_Lymph_Node.h5seurat"), Results)
}

sc.tab.ln <- SeuratDisk::LoadH5Seurat(paste0(Results, "TS_Lymph_Node.h5seurat"),
                                      assays = "RNA")


# leaving out the samples generated using 'smartseq2' protocol
sc.tab.ln <- subset(x = sc.tab.ln, subset = method == "10X")


# Seurat expects nCount_RNA and nFeature_RNA columns in metadata referring to
# total UMI detected and total unique genes detected in a cell respectively
# Changing the names to adher to Seurat standards
sc.tab.ln@meta.data$nCount_RNA <- sc.tab.ln@meta.data$n_counts_UMIs
sc.tab.ln@meta.data$nFeature_RNA <- sc.tab.ln@meta.data$n_genes

sc.tab.ln@meta.data$n_counts_UMIs <- NULL
sc.tab.ln@meta.data$n_genes <- NULL

# Add number of genes per UMI for each cell to metadata
sc.tab.ln$log10GenesPerUMI <-
  log10(sc.tab.ln$nFeature_RNA) / log10(sc.tab.ln$nCount_RNA)

# Compute percent mito ratio
sc.tab.ln$mitoRatio <- PercentageFeatureSet(object = sc.tab.ln, pattern = "^MT-", assay = "RNA")
sc.tab.ln$mitoRatio <- sc.tab.ln@meta.data$mitoRatio / 100


# annotate entire tabula sapiens lymph node dataset using SingleR package
sc.tab.ln0 <- annotate_clusters(sc.tab.ln)

tmp0 <- sc.tab.ln0[, !is.na(sc.tab.ln0@meta.data[, "blue.main"])]

png(file = paste0(Results, "suppl_fig_1a.png"), res = 450,
    width = 7.2, height = 4.8, units = "in")
DimPlot(tmp0,
        reduction = "umap",
        cols = get_pallete(length(unique(tmp0$blue.main))),
        label = T,
        repel = T,
        label.size = 2,
        label.box = T,
        group.by = "blue.main") +
  ggtitle("Tabula sapiens lymph node full dataset") +
  NoLegend()
dev.off()
rm(tmp0, sc.tab.ln0)


#### splitting tabula sapiens dataset ####

# Splitting tabula sapiens data donor wise
Idents(sc.tab.ln) <- "donor"

sc.tab1 <- subset(sc.tab.ln, idents = "TSP2")
sc.tab2 <- subset(sc.tab.ln, idents = "TSP7")
sc.tab3 <- subset(sc.tab.ln, idents = "TSP14")


# subsetting the tabula sapiens lymph node dataset based on 10X's profiling
cell_names <- colnames(sc.tab3)

prime5 <- list() # cells with 10X's 5-prime profiling
prime3 <- list() # cells with 10X's 3-prime profiling

# Split the list based on the 10X's profiling
for (value in cell_names) {
  if (grepl("5Prime", value)) {
    prime5 <- c(prime5, value)
  } else {
    prime3 <- c(prime3, value)
  }
}

sc.5prime <- sc.tab3[,colnames(sc.tab3) %in% unlist(prime5)]
sc.3prime <- sc.tab3[,colnames(sc.tab3) %in% unlist(prime3)]


#### Filtering tabula sapiens data ####
# cell-level filtering
sc.tab.ln.fil1 <- subset(
  x = sc.tab1,
  subset = (nCount_RNA >= 2300) &
    (nFeature_RNA >= 750) &
    (log10GenesPerUMI > 0.80) &
    (mitoRatio < 0.30))

# Process tabula sapiens data
sc.tab.ln.nor1 <- process_data(sc.tab.ln.fil1)

# annotate tabula sapiens data with SingleR annotation package
sc.tab.ln.nor1 <- annotate_clusters(sc.tab.ln.nor1)


# cell-level filtering
sc.tab.ln.fil2 <- subset(
  x = sc.tab2,
  subset = (nCount_RNA >= 2300) &
    (nFeature_RNA >= 750) &
    (log10GenesPerUMI > 0.80) &
    (mitoRatio < 0.30))

# Process tabula sapiens data
sc.tab.ln.nor2 <- process_data(sc.tab.ln.fil2)

# annotate tabula sapiens data with SingleR annotation package
sc.tab.ln.nor2 <- annotate_clusters(sc.tab.ln.nor2)


# cell-level filtering
sc.tab.ln.fil3 <- subset(
  x = sc.5prime,
  subset = (nCount_RNA >= 2300) &
    (nFeature_RNA >= 750) &
    (log10GenesPerUMI > 0.80) &
    (mitoRatio < 0.30))

# Process tabula sapiens data
sc.tab.ln.nor3 <- process_data(sc.tab.ln.fil3)

# annotate tabula sapiens data with SingleR annotation package
sc.tab.ln.nor3 <- annotate_clusters(sc.tab.ln.nor3)

# cell-level filtering
sc.tab.ln.fil4 <- subset(
  x = sc.3prime,
  subset = (nCount_RNA >= 2300) &
    (nFeature_RNA >= 750) &
    (log10GenesPerUMI > 0.80) &
    (mitoRatio < 0.30))

# Process tabula sapiens data
sc.tab.ln.nor4 <- process_data(sc.tab.ln.fil4)

# annotate tabula sapiens data with SingleR annotation package
sc.tab.ln.nor4 <- annotate_clusters(sc.tab.ln.nor4)


#### In-house LNSC dataset ####
plus <- Read10X(data.dir = paste0(Data, "in_house_lnsc/plusDataset"))
plus <- CreateSeuratObject(counts = plus, project = "plus")

minus <- Read10X(data.dir = paste0(Data, "in_house_lnsc/minusDataset"))
minus <- CreateSeuratObject(counts = minus, project = "minus")

sc.in.house <- merge(plus,
                     y = minus,
                     add.cell.ids = c("plus", "minus"),
                     project = "combined")
rm(plus, minus)

# Add number of genes per UMI for each cell to metadata
sc.in.house$log10GenesPerUMI <-
  log10(sc.in.house$nFeature_RNA) / log10(sc.in.house$nCount_RNA)

# Compute percent mito ratio
sc.in.house$mitoRatio <- PercentageFeatureSet(object = sc.in.house,
                                              pattern = "^MT-",
                                              assay = "RNA")
sc.in.house$mitoRatio <- sc.in.house@meta.data$mitoRatio / 100
range(sc.in.house$mitoRatio)


#### Filtering in-house data ####
# cell-level filtering
sc.in.house.fil <- subset(x = sc.in.house,
                          subset = (nCount_RNA >= 500) &
                            (nFeature_RNA >= 250) &
                            (log10GenesPerUMI > 0.70) &
                            (mitoRatio < 0.20))


# Process in-house lymph node data
sc.in.house.nor <- process_data(sc.in.house.fil)

# annotate in-house LNSC using SingleR package
sc.in.house.nor <- annotate_clusters(sc.in.house.nor)

# UMAP representation of the in-house lymph node full dataset without NA
tmp <- sc.in.house.nor[, !is.na(sc.in.house.nor@meta.data[, "blue.main"])]

png(file = paste0(Results, "suppl_fig_1b.png"), res = 450,
    width = 7.2, height = 4.8, units = "in")
p1 <- DimPlot(tmp,
              reduction = "umap",
              cols = get_pallete(length(unique(tmp$blue.main))),
              label = T,
              repel = T,
              label.size = 2,
              label.box = T,
              group.by = "blue.main") +
  ggtitle("Lymph node stromal cells full dataset") +
  NoLegend()
print(p1)
dev.off()
rm(tmp)


#### Integrated data ####
# combining both datasets
sc.datas <- list(sc.in.house.nor, sc.tab.ln.nor1, sc.tab.ln.nor2,
                 sc.tab.ln.nor3, sc.tab.ln.nor4)

# Select the most variable features to use for integration
features <- SelectIntegrationFeatures(object.list = sc.datas, nfeatures = 3000)

# Find best buddies - can take a while to run
sc.anchors <- FindIntegrationAnchors(object.list = sc.datas,
                                     anchor.features = features)


# Integrate across conditions
sc.combined <- IntegrateData(anchorset = sc.anchors)


sc.combined <- ScaleData(sc.combined, verbose = F) %>%
  RunPCA(verbose = F) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = F)

# Annotate integrated data
sc.combined <- annotate_clusters(sc.combined)

saveRDS(sc.combined, file = paste0(Results, "sc.combined.nor.rds"))


#### Downsample on features ####
counts <- sc.combined@assays$RNA@counts

# percent value for each cell where the gene is expressed
gene.percent.exp <- rowMeans(counts > 0) * 100

# to select genes detected in at least 3% of cells
genes.filter.1 <- names(gene.percent.exp[gene.percent.exp > 3])

# to select genes detected in few cells between 0.01% - 3%
genes.filter.2 <- names(gene.percent.exp[gene.percent.exp < 3 &
                                           gene.percent.exp > 0.01])

# mean expression of genes across non-zero cells
nonz.mean <- (rowSums(counts) / rowSums(counts > 0))

# to select genes with mean expression more than 2 across non-zero cells
genes.filter.3 <- names(nonz.mean[nonz.mean > 2])

genes.filter <- union(genes.filter.1, intersect(genes.filter.2, genes.filter.3))


# count matrix with filtered genes on the integrated assay
sc.data <-
  CreateSeuratObject(counts = sc.combined@assays$RNA@counts[genes.filter, ],
                     meta.data = sc.combined@meta.data)


sc.data <- FindVariableFeatures(sc.data,
                                selection.method = "vst",
                                nfeatures = 3000,
                                verbose = F) %>%
  ScaleData(verbose = F) %>%
  RunPCA(npcs = 30, verbose = F) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = F)


# removing cells with NA annotations
sc.data.subset <- sc.data[, !is.na(sc.data@meta.data[, "blue.main"])]

png(paste0(Results, "suppl_fig_1c.png"),
    res = 450, width = 8.4, height = 6, units = "in")
p2 <- DimPlot(sc.data.subset,
              reduction = "umap",
              cols = get_pallete(length(unique(sc.data.subset$blue.main))),
              label = T,
              repel = T,
              label.size = 2.5,
              label.box = T,
              group.by = "blue.main") +
  ggtitle("Integrated complete scRNA-seq dataset") +
  NoLegend()
print(p2)
dev.off()


#### Downsample on number of cells ####

# choosing idents with cells equal or more than Neutrophils
nr_cells_cutoff <- unname(sort(table(sc.data.subset$blue.main))["Neutrophils"])
idents.use <- data.frame(sort(table(sc.data.subset$blue.main))) %>%
  filter(Freq >= nr_cells_cutoff) %>%
  dplyr::select(Var1)

Idents(sc.data.subset) <- "blue.main"
sc.data.pro <- subset(sc.data.subset, idents = idents.use$Var1)

# downsampling cells per cell type (cluster)
# WhichCells uses sample() function internally to select random cells
cell.lists <- WhichCells(sc.data.pro, downsample = 300, seed = 15)
length(cell.lists)

sc.data.pro <- sc.data.pro[, cell.lists]


#### updating cell types names suitable for programming standards
current.cluster.ids <- c("Adipocytes",
                         "B-cells",
                         "CD4+ T-cells",
                         "CD8+ T-cells",
                         "Endothelial cells",
                         "Fibroblasts",
                         "HSC",
                         "Macrophages",
                         "Monocytes",
                         "Myocytes",
                         "Neutrophils",
                         "NK cells",
                         "Skeletal muscle")
new.cluster.ids <- c("Adipocytes",
                     "B_cells",
                     "CD4_T_cells",
                     "CD8_T_cells",
                     "Endothelial_cells",
                     "Fibroblasts",
                     "HSC",
                     "Macrophages",
                     "Monocytes",
                     "Myocytes",
                     "Neutrophils",
                     "NK_cells",
                     "Skeletal_muscle")
sc.data.pro@meta.data$blue.main <- plyr::mapvalues(x = sc.data.pro@meta.data$blue.main,
                                                   from = current.cluster.ids,
                                                   to = new.cluster.ids)


# to get the exact same colors for cell types as that of full dataset
color_pal <- c("#1B9E77", "#7570B3", "#E7298A", "#66A61E", "#666666", "#B2DF8A",
               "#33A02C", "#FB9A99", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99",
               "#B15928")

png(paste0(Results, "suppl_fig_1d.png"),
    res = 450, width = 8.4, height = 6, units = "in")
p3 <- DimPlot(sc.data.pro,
              reduction = "umap",
              cols = color_pal,
              label = T,
              repel = T,
              label.size = 2.5,
              label.box = T,
              group.by = "blue.main") +
  ggtitle("Integrated downsampled scRNA-seq reference dataset") +
  NoLegend()
print(p3)
dev.off()
rm(color_pal)


saveRDS(sc.data.pro, file = paste0(Results, "sc.ref.data.rds"))

