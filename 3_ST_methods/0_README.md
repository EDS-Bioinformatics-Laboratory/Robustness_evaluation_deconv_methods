# Reference-based cell type deconvolution methods



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


**Software environment:** 

- R version : 4.1.2 (2021-11-01)
- Rstudio IDE: RStudio 2023.03.0+386 "Cherry Blossom"



####**Conceptual description of methodology**


In our systemic analysis, we have included reference-based cell type deconvolution methods tailored for spatial transcriptomics data as well as bulk RNA-seq data.

* **CARD**: This method uses a non-negative matrix factorization model to use the cell type-specific gene expressions from scRNA-seq reference data. The method accommodates the spatial correlation structure in cell type composition across tissue by conditional autoregressive (CAR) modeling. Thus, we have used the spatial coordinates for the simulated ST data as one of the inputs. For the analysis, we follow the guidelines from the [tutorial](https://yingma0107.github.io/CARD/documentation/04_CARD_Example.html) provided by the method's developers. All the parameters are set to default as suggested in the tutorial.

* **Cell2Location**: This method uses a Bayesian model to estimate the absolute abundance of cell types at each spot by decomposing the spatial expression count matrix into a predefined set of reference cell type signatures. For the analysis, we follow the scvi-tools implementation [guidelines](https://docs.scvi-tools.org/en/0.15.1/tutorials/notebooks/cell2location_lymph_node_spatial_tutorial.html) of the Cell2Location model. We filtered the genes again with the parameter value set to cell\_count\_cutoff = 5, cell\_percentage\_cutoff = 0.03, nonz\_mean\_cutoff = 1.12. The latter two cutoffs have been applied to scRNA-seq reference data earlier; the result remains the same. For training the reference scRNA-seq data, we used 1000 epochs and for the rest of the parameters default values were used. For creating the model for spatial data, we set the N\_cells\_per\_location = 13 based on the average number of cells present per spot in our simulated ST data; all other parameter values were used as suggested by the developers. While for training the ST data model we used 10000 epochs compare to 30000 as suggested by the developers. We used the 5% quantile of the posterior distribution as proposed in the guidelines that represented the cell type abundance in each spot. The cell type abundance values for each spot were normalized further to use in the downstream comparative analysis.

* **RCTD**: This method uses a supervised learning approach to decompose RNA sequencing mixtures into single cell types and enables the assignment of the cell types in a spot. The RCTD method utilizes the spatial coordinates of the simulated ST data. We have followed the [guidelines](https://raw.githack.com/dmcable/spacexr/master/vignettes/spatial-transcriptomics.html) provided by the developers for applying RCTD to spatial transcriptomics data. For creating the reference data constructor and spatial data constructor, we stick with the default parameter values. While creating and executing the RCTD object, we used the default parameters, and the number of cores was set to 8 for parallel processing. The RCTD was used in the full mode (i.e., doublet_mode = ‘full’) where the model can fit any number of cell types on each spot.

* **Seurat**: We have used the Seurat V4 workflow in our study. It uses the weighted-nearest neighbor algorithm which is based on the unsupervised strategy of learning cell-specific modality weights. The two modalities i.e., scRNA-seq and spatial transcriptomics are integrated to find the cell type proportions present in each spot using the [guidelines](https://satijalab.org/seurat/articles/integration_introduction.html) provided by the method developers. We have used the standard normalization procedure with 3000 highly variable features as demonstrated in Seurat V3 to normalize the scRNA-seq and ST data before integration. To integrate we first find the anchors between a scRNA-seq and ST data using ‘FindTransferAnchors’ functionality which utilizes the PCA space constructed using the scRNA-seq data. Later we used these anchors to detect the cell types present in each spot by transferring the data from the scRNA-seq modality to the ST modality using ‘TransferData’ functionality. All functions are used with default parameter values unless stated otherwise.

* **SPOTlight**: This method is based on non-negative matrix factorization (NMF) regression. It also utilizes the non-negative least squares (NNLS) to populate the coefficient matrices with cell type marker genes. The coefficient and basis matrices are initiated by adding prior information to the model. We have followed the guidelines provided in the [vignette](https://github.com/MarcElosua/SPOTlight/blob/main/vignettes/SPOTlight_kidney.Rmd) by the method's developer. We used all the cells in scRNA-seq reference data to get the best of the method. The earlier calculated 3000 highly variable genes in our integrated scRNA-seq reference data were used without any change. We changed the marker gene scores cut off from 0.8 to 0.5 (mean.AUC > 0.5) to ensure more relevant marker genes suitable to our reference data. The SPOTlight deconvolution algorithm has been executed with the default parameters for the rest of the analysis.

* **Stereoscope**: This method uses a model-based approach, where a negative binomial model is optimized for a gene and a cell type-specific parameter. For the analysis, we follow the scvi-tools implementation [guidelines](https://docs.scvi-tools.org/en/stable/tutorials/notebooks/stereoscope_heart_LV_tutorial.html) of the Stereoscope model. We skip the pre-processing of the scRNA-seq data since it has already been performed. Further, we only select 5000 highly variable genes instead of 7000 as mentioned in the documentation. The rest of the hyper-parameters in the estimation of expression signature from scRNA-seq data and inferred proportion of ST data has been set to default.

* **MuSiC**: This method was developed to estimate the cell type proportions in bulk RNA-seq data with the help of weighted non-negative least squares (W-NNLS). The method uses cross-subject and cross-cell consistencies to guard against bias in subject selection and cell capture in scRNA-seq data. We have followed the guidelines provided in the [vignette tutorial](https://xuranw.github.io/MuSiC/articles/MuSiC.html) by the method's developers. We estimated the cell type proportions without pre-grouping on cell types. As suggested all the genes are used during the execution of the method without any prior selection of marker genes; this ensures the common genes are used for the cell type estimation. All the other parameter values are set to default values. We assume each spot in the ST data as one sample, resulting in 1600 samples of bulk RNA-seq data.

* **SCDC**: This method adopts an ENSEMBLE approach to integrate deconvolution results from across datasets. The method supports multiple scRNA-seq reference datasets by addressing the problem of batch effect. We have used the ‘SCDC\_prop\_ONE’ functionality which is suitable for our single scRNA-seq reference data. SCDC method utilizes the W-NNLS framework proposed by [9] with a slightly different approach to calculating the basis matrix. We have followed the guidelines suggested by the developers in the [vignette tutorial](https://meichendong.github.io/SCDC/articles/SCDC.html). We specifically followed the three-cell-line mixture data analysis tutorial. Here also we assume each spot in the ST data as one sample, resulting in 1600 samples of what we call bulk RNA-seq data. All the parameters were used with the default values.


