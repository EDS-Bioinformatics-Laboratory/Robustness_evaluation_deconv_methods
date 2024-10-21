# Post analysis work and manuscript figures



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

- R version : 4.1.2 (2021-11-01)
- Rstudio IDE: RStudio 2023.03.0+386 "Cherry Blossom"



**Conceptual description of methodology**

We used JSD and RMSE metrics to analyse the reference-based cell type deconvolution method results, and here we have R scripts to generate the most significant results in the project which are also included in the manuscript.


**Is a script and data available to regenerate the important figures?**

The standalone R scripts generates resulting data and figures during the execution.



