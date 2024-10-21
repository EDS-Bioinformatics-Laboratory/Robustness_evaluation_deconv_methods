# Generating simulated spatial transcriptomics datasets



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


We generate a single spot and apply the same algorithm to the rest of the spots. In our simulated spatial transcriptomics, we consider each spot an independent entity inferring no biological relevance among neighboring spots. This assumption suits our goal of establishing a ground truth and only focuses on the deconvolution of a spot for cell type proportion using the reference scRNA-seq data.

* Decide the number of cells present in a spot randomly between a range

* Decide the number of cell types present in that spot. This random selection of numbers follows uniform distribution when we look at all the spots.

* Calculate the probability of each cell type existing in that spot based on the number of cell types and cells.

* Using the probability of each cell type in a spot, we calculate the number of cells per cell type with the help of multinomial distribution. This distribution of proportions of cell types will serve as the *ground truth* in the analysis.

* The mRNA count for each spot is calculated using the following equation,		
$$
mRNA\ count\ for\ spot\ i = \sum_{k = 1}^{N}{{R_{avg, k}} * W_{i, k}}
$$
	here, N is the number of cell types in the spot, R<sub>avg, k</sub> is the mRNA count for cell type k, and W<sub>i, k</sub> is the number of cells for cell type k in the spot. Ravg is calculated for each spot using five randomly selected cells of a cell type. Choosing five random cells of a cell type and considering their average for each spot cancels the possibility of selecting the same cell over and over and generating a biased simulated dataset. In the scenarios where mRNA counts exceed a limit of 30000 for any spot, it is down-sampled to a random number between 25000-30000 to harmonize with the spot of actual 10x Genomics visium-like spatial transcriptomics data.

* We generate dummy spatial coordinates for the deconvolution methods, which rely on the coordinates of spots in the tissue slide.



