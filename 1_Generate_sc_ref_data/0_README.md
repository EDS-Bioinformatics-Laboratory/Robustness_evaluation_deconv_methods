# Generate single-cell reference dataset



*<u>**Instructions**</u>*

* *Remove instructions once you finalized the content of this file.*
* *Subdirectories*
  * ***Code**. Contains the (in-house developed) software for the computational analysis.*
  * ***CodeDocumentation**. Any documentation that is available about the software such as requirement specifications, software design, source code documentation, testing requirements, and end-user instructions. Alternatively, part of this documentation can be placed in the code itself, or in the code 0_README.md file.*
  * ***Data**. Raw, meta, and pre-processed data.*
  * ***Notebooks**. Contains interactive notebooks (e.g, R or Jupyter notebooks), which may include data and results. Minimize the size of the Notebook prior to pushing to GitHub.*
  * ***Results**. (Intermediate) results such as figures and tables from the computational analysis. Additional sub-directories within /Results, to organize the results, are allowed.* 
  * ***Settings**. Contains any file with input parameters for e.g., statistical analysis and computational simulations. If only few parameters or parameter files are required then these might also be placed in the \Code directory or directly in the code itself.*



*<u>Manual steps:</u>*

***Avoid** manual data manipulation steps. Instead of manually changing data (e.g., format conversion, filtering), aim to make use of small custom scripts. Manual steps are a main cause for irreproducible results.*  

* *Having an executable description of the sequence of steps taken helps to reproduce results (e.g., Snakemake (https://snakemake.github.io/) or pytask (https://github.com/pytask-dev/pytask)).*



*<u>Intermediate results and hierarchical output</u>*

*Record intermediate results (preferable in a standardized format). Generate hierarchical analysis output, allowing layers of increasing detail to be inspected. This can reveal discrepancies toward what is assumed, and can in this way uncover bugs or faulty interpretation that are not apparent in the final results. It also allows any inconsistency to be tracked to the step where the problem occurs.  It also allows critical examination of the full process behind a result.*



*<u>Regeneration of figures and tables:</u>* 

*For any figure or table that ends up in a publication, report, or presentation at meeting, the underlying data and a stand-alone piece of code should be available to regenerate the figure. It also allows easy modification of a figure and to retrieve the data of the figure (instead of having to redo a complete analysis). Equally important, the data of the figure can be further analysed or inspected.*



*== END INSTRUCTIONS ==*



**Operating System(s) / version(s) used during development (and testing):**

MacOS Ventura 13.4.1

**Specific hardware requirements:** NA

**Software environment:** 

* R version : 4.1.2 (2021-11-01)
* Rstudio IDE: RStudio 2023.03.0+386 "Cherry Blossom"



####**Conceptual description of methodology**


* We used two independent scRNA-seq datasets in the project to create a comprehensive scRNA-seq reference dataset using an integratio strategy.

* The integration strategy used during the analysis uses canonical correlation analysis and mutual nearest neighbor approaches to identify the cells in shared space across the dataset, which serves as anchors to guide the dataset integration.
* Both the scRNA-seq datasets are processed differently until they are integrated. Below the details of data filtering, reduction, normalization.

	* 	**scRNA-seq lymph node data from tabula sapiens consortium**
		*  	Filtering: Joint effects of the below metrics 
			*   minimum 2300 molecules within a cell (nCount_RNA > 2300)
			*   minimum 750 genes per cell (nFeature_RNA > 750)
			*   less than 20% mitochondrial genes per cell
			*   genes detected per UMI > 0.80 (on logarithmic scale) gives idea about the complexity of the data (higher the values, complex the data). Low complexity refers to the contamination
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
			*   genes detected per UMI > 0.70 (on logarithmic scale) gives idea about the complexity of the data (higher the values, complex the data). Low complexity refers to the contamination
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

* For the integration of these two distinct datasets we followed the instructions in the [integration vignette](https://satijalab.org/seurat/articles/integration_introduction.html) suggested by the developers with a little deviation in number of integration anchors as 3000 (default = 2000).

* The integrated reference data is further scaled and reduced dimensions were calculated with settings similar to the either of scRNA-seq dataset.

* For the integrated scRNA-seq data filtering:
	* for *filtering on features*, we utilized genes expressed in at least 3% of the cells and genes expressed between 0.01% - 3%, but with the mean expression of genes across non-zero cells greater than 2.
	* for *filtering on cells*, we annotated the cells using functionality offered by SingleR (R package - version 1.8.1) and selected celltypes with minimum 100 cells and to maintain a homogeneous mixture of cells at max 300 cells per celltypes.
