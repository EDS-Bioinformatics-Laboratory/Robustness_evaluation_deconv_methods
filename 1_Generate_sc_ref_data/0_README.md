# Generate a single-cell reference dataset


**Operating System(s) / version(s) used during development (and testing):**

MacOS Ventura 13.4.1

**Specific hardware requirements:** NA

**Software environment:** 

* R version : 4.1.2 (2021-11-01)
* Rstudio IDE: RStudio 2023.03.0+386 "Cherry Blossom"



####**Conceptual description of methodology**


* In the project, we used two independent scRNA-seq datasets (Tabula sapiens lymph node dataset and lymph node stromal cells dataset GEO accession ID: #######) to create a comprehensive scRNA-seq reference dataset using an integration strategy.

* The integration strategy used during the analysis uses canonical correlation analysis and mutual nearest neighbour approaches to identify the cells in shared space across the dataset. These serve as anchors to guide the dataset integration.

* Both the scRNA-seq datasets are processed differently until integrated below the details of data filtering, reduction, and normalisation.

	* 	**scRNA-seq lymph node data from tabula sapiens consortium**: The tabula sapiens lymph node data comprises three different donors, so we split it into three independent datasets; further, one of the donor datasets is split into two platform-wise. Later, we applied the same filtering, reduction and normalisation strategy.
	
		*  	Filtering: Joint effects of the below metrics 
			*   minimum 2300 molecules within a cell (nCount_RNA > 2300)
			*   minimum 750 genes per cell (nFeature_RNA > 750)
			*   less than 20% mitochondrial genes per cell
			*   genes detected per UMI > 0.80 (on a logarithmic scale) give an idea about the complexity of the data (the higher the values, the more complex the data). Low complexity refers to the contamination
		*   Normalization: Normalize the count data
			*   method = LogNormalize
			*   scale.factor = 10000
		*	Variable Features: Identifies features that are outliers on a 'mean variability plot'
			*	selection method = "vst"
			*	number of features = 3000
		*	Scaling: Scales and centers features in the dataset
			*	default settings offered by Seurat workflow
		*	Dimensionality reduction: Using PCA and UMAP
			* number of principal components = 30

	* 	**scRNA-seq lymph node stromal cell data**
		*  	Filtering: Joint effects of the below metrics 
			*   minimum 500 molecules within a cell (nCount_RNA > 500)
			*   minimum 250 genes per cell (nFeature_RNA > 250)
			*   less than 20% mitochondrial genes per cell
			*   genes detected per UMI > 0.70 (on a logarithmic scale) give an idea about the complexity of the data (the higher the values, the more complex the data). Low complexity refers to the contamination
		*   Normalization: Normalize the count data
			*   method = LogNormalize
			*   scale.factor = 10000
		*	Variable Features: Identifies features that are outliers on a 'mean variability plot'
			*	selection method = "vst"
			*	number of features = 3000
		*	Scaling: Scales and centers features in the dataset
			*	default settings offered by Seurat workflow
		*	Dimensionality reduction: Using PCA and UMAP
			* number of principal components = 30

* For the integration of these two distinct datasets, we followed the instructions in the [integration vignette](https://satijalab.org/seurat/articles/integration_introduction.html) suggested by the developers, with a slight deviation in the number of integration anchors to 3000 (default = 2000).

* The integrated reference data is further scaled, and reduced dimensions were calculated with settings similar to those of the scRNA-seq dataset.

* For the integrated scRNA-seq data filtering:
	* for *filtering on features*, we utilized genes expressed in at least 3% of the cells and genes expressed between 0.01% - 3%, but with the mean expression of genes across non-zero cells greater than 2.
	* for *filtering on cells*, we annotated the cells using functionality offered by SingleR (R package - V1.8.1) and selected cell types with a minimum of 100 cells and maintained a homogeneous mixture of cells at a maximum of 300 cells per cell types.
