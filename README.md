# Systematic evaluation of robustness to cell type mismatch of deconvolution methods for spatial transcriptomics data

## Table of contents

* [Overview](#overview)
* [Contents of the current directory](#contents-of-the-current-directory)
* [Installation steps](#installation-steps)
* [Reproducing manuscript results](#how-can-you-reproduce-the-results-in-the-manuscript)
	- [Module 1: Generate single-cell reference datasets](#module-1-generate-single-cell-reference-datasets)
	- [Module 2: Generate simulated ST data](#module-2-generate-simulated-st-data)
 	- [Module 3: Deconvolution method results](#module-3-deconvolution-method-results)
  	- [Module 4: Analyse deconvolution results](#module-4-analyse-deconvolution-results)
 	- [Module 5: Post analysis results](#module-5-post-analysis-results)	 	
* [Known issues](#known-issues)
* [Appendix A](#appendix-a)
<br>

## Overview
&emsp;&emsp;&emsp;&emsp; Spatial transcriptomics approaches based on sequencing (barcode-based, e.g., 10x Visium) preserve spatial information but with limited cellular resolution. On the other hand, single-cell RNA-sequencing (scRNA-seq) techniques provide single-cell resolution but lose spatial resolution because of the tissue dissociation step during the scRNA-seq experimental procedure. With these complementary strengths in mind, computational methods have been developed to combine scRNA-seq and spatial transcriptomics data. These approaches use deconvolution to identify cell types and their respective proportions at each location in spatial transcriptomics data with the aid of a scRNA-seq reference dataset. Some suggest that deconvolution approaches are sensitive to the absence of cell type(s) in the single-cell reference dataset, a problem referred to as *cell type mismatch*. Here, we systematically evaluated the robustness of deconvolution methods to cell type mismatch tailored for spatial transcriptomics data.
<br>

## Contents of the current directory

Note that in the complete sFSS, this directory corresponds to the Processing directory. 

* **0_SoftwareEnvironment**:
This directory enlists the software environment specifications used for various programming languages and/or platforms during the project.

* **1\_Generate\_sc\_ref\_data**: The directory comprises scripts, results and settings to generate the integrated scRNA-seq dataset from 2 distinct and complementary scRNA-seq datasets. The final result **`sc.ref.data.rds`** is used as a basis to generate various reference dataset versions based on cell type removal scenarios and to generate simulated spatial transcriptomics datasets.

* **2\_Simulating\_ST\_data**: Comprised of scripts, results and settings to generate the sequencing-based simulated spatial transcriptomics datasets varying in number of cells & cell types present per spatial location (spot). We created three simulated ST datasets using the algorithm developed for the analysis. <br>
&emsp;&emsp;&emsp; a. **ST1**: 4-8 cell types and 10-15 cells per spot; <br>
&emsp;&emsp;&emsp; b. **ST2**: 1-5 cell types and 10-15 cells per spot; <br>
&emsp;&emsp;&emsp; c. **ST3**: 1-5 cell types and 3-7 cells per spot.

* **3\_ST\_methods**: Comprised of standalone R/Python scripts, one for each deconvolution method and shell/batch scripts to execute the R or Python scripts in parallel for multiple instances in a removal scenario, provided the required computational power is available.
	- **R-based**: CARD, RCTD, SCDC, MuSiC, Seurat, SPOTlight
	- **Python-based**: cell2location, Stereoscope (GPU recommended) <br>
 
	&emsp;&emsp;&emsp;The deconvolution results for each removal scenario are available in the respective directory (for instance, the _**rm1**_ directory refers to the removal scenario for removing one cell type from the scrna-seq reference data). See below the overview of removal scenarios and total reference datasets in each scenario. <br>
	- **rm0** - removal of _**no**_ cell types from reference data (1 dataset) <br>
	- **rm1** - removal of _**one**_ cell type from reference data (13 datasets) <br>
	- **rm2** - removal of _**two**_ cell types from reference data (5 datasets) <br>
	- **rm3** - removal of _**three**_ cell types from reference data (5 datasets) <br>
	- **rm5** - removal of _**five**_ cell types from reference data (1 dataset) <br>
	- **rm10** - removal of _**ten**_ cell types from reference data (5 datasets) <br>
	- **rm11** - removal of _**eleven**_ cell types from reference data (5 datasets)

* **4\_Analysis\_results**: Comprised of scripts to calculate similarity metrics like JSD and RMSE to understand the performance of deconvolution methods for various cell type removal scenarios compared to a baseline with no cell type missing. Cell type reassignment metrics calculate the assignment of missing cell type proportions from the reference dataset.

* **5\_Post\_analysis\_work**: The directory comprises scripts to generate the final resulting plots in the project; a few plots are included in the manuscript's main text, while others were included in the supplementary information.

* **Data** (only available in sFSS, not on GitHub): Comprised of pre-processed data downloaded from a public data-sharing platform. One of the datasets is procured from the group of Lisa van Baarsen and is available on GEO as well with GEO accession ID - ####. 
* **renv**: The project uses the renv functionality and the directory comprises the essential files (activate.R, settings.json) and directories for the renv setup.

* **Renv_setup.R**: This script sets up the R environment for the project work. Details about executing are available in the file as a header and in the comments. Also included under the 'How to reproduce the results' section further down in this document.<br>

* **renv.lock**: R environment lock (metadata) file comprising package details. This file will be used by the *Renv_setup.R* script to download and install the correct package and its version to recreate the computing environment.

* **.Rprofile**: Essential R profile script to set up the correct R profile while using the renv infrastructure.

* **environment.yml**: A metadata file comprising Python packages installed in the virtual conda environment.
<br>

## Installation steps

Running the code requires a Linux OS and the R, RTools and Python versions indicated below. 

#### Linux / Ubuntu distribution

1. Install Anaconda ([easy guide for installation](https://phoenixnap.com/kb/install-anaconda-ubuntu))
2. Install R version 4.1.2 ([download](https://cran.rstudio.com/src/base/R-4/R-4.1.2.tar.gz)). Instructions on [how to install](https://docs.posit.co/resources/install-r-source.html)
3. Set up conda environment with Python version 3.9.7 using the command: <t>`conda env create -f environment.yml`
4. Set up renv setup details
	- Execute `Renv_setup.R` script from where it resides currently to initialise the renv infrastructure for the R project (ignore the warning messages) using `Rscript Renv_setup.R arg1` command. <br>
	Expected command line argument (*arg1*) is **(classic) GITHUB_PAT token** ([more details](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens)).<br>
	*The required files for renv setup (renv.lock, .Rprofile, renv/activate.R, renv/settings.json) should already be in your R project directory; If not, ensure it resides in the same directory as the Renv_setup.R file.*
> [!NOTE]
> Due to the platform (os) dependencies, a collaborator can still experience minor discrepancies in the results while using the conda/renv functionality.

Installation steps for Windows OS can be found in [Appendix A](./#appendix-a). Note however that using Windows only R-based deconvolution methods can be executed and that therefore the results of Python-based methods (cell2location, SPOTlight) cannot be reproduced. 
<br>
<br>

## How can you reproduce the results in the manuscript?
> [!WARNING] 
> The following instructions are for reproducing the results using the command-line terminal/prompt, not RStudio or Jupyter Notebook.

### Module 1: Generate single-cell reference datasets
#### Navigate to "1\_Generate\_sc\_ref_data/Code/" directory
* Generating a single-cell reference dataset from multiple single-cell datasets, UMAP representations included in the supplementary information of the manuscript using the commands as below;
	
	&emsp;&emsp;&emsp;&emsp;`Rscript Init_env.R` <br>
	&emsp;&emsp;&emsp;&emsp;`Rscript SC_ref_data.R` 
	
	> Note: this script takes a significant amount of time (~2-3 hours on a machine with 16GB RAM) to execute. Ignore the warnings during the execution
	
* Generating single-cell reference datasets for all removal scenarios based on cell type removal using the commands as below;
	
	&emsp;&emsp;&emsp;&emsp;`Rscript SC_ref_data_scenarios.R`
	> Note: this script requires the single-cell reference dataset generated in the previous step.

### Module 2: Generate simulated ST data
#### Navigate to "2\_Simulating\_ST\_data/Code/" directory
* Generating simulated spatial transcriptomics datasets using single-cell reference data, the three datasets vary by the number of cells and cell types present per spatial location using the commands as below;

	&emsp;&emsp;&emsp;&emsp;`Rscript Init_env.R` <br>
	&emsp;&emsp;&emsp;&emsp;`Rscript Generate_ST_data.R` *arg1 arg2 arg3 arg4 arg5*	
	> &emsp;&emsp;&emsp;Note: `Generate_ST_data.R` script expects five command line arguments in the order mentioned below;<br>
	>	&emsp;&emsp;&emsp;&emsp; arg1 - min number of cell types per spot <br>
	>	&emsp;&emsp;&emsp;&emsp; arg2 - max number of cell types per spot <br>
	>	&emsp;&emsp;&emsp;&emsp; arg3 - min number of cells per spot <br>
	>	&emsp;&emsp;&emsp;&emsp; arg4 - max number of cells per spot <br>
	>	&emsp;&emsp;&emsp;&emsp; arg5 - index of simulated ST dataset [*options: 1, 2, 3* ] <br>
	
		To reproduce the results in the manuscript, the following command line input generates the first simulated ST dataset:
		Rscript Init_env.R
		Rscript Generate_ST_data.R 4 8 10 15 1

### Module 3: Deconvolution method results
#### Navigate to "3\_ST\_methods/Code/" directory 
* The shell scripts execute all deconvolution methods to predict cell type proportions simultaneously for all removal scenarios and multiple reference datasets for each scenario. <br>

	Shell scripts to execute R and Python-based deconvolution methods: <br>
	&emsp;&emsp;&emsp;&emsp;`Rscript Init_env.R` <br>
	&emsp;&emsp;&emsp;&emsp;`Execute_R_based_methods.sh` *arg1 arg2 arg3*, <br>
	&emsp;&emsp;&emsp;&emsp;`Execute_python_based_methods.sh` *arg1 arg2 arg3*, <br>

	
	Please make sure to execute only one script at a time and wait for the results. For every removal scenario, each single-cell reference will have one result (for instance, rm1 will have 13 results for each method) for every method.
	
	 *If your machine crashes, please adapt the shell script accordingly (for instance, start fewer jobs in parallel). Also, verify the number of results generated before moving to the next scenario. If a combination of SC & ST datasets fails for a method, please run the standalone script for that pair.* <br>

	> Note: The script to execute R-based/Python-based methods executes each method in parallel for a provided number of reference datasets and removal scenarios; the required command line arguments are as below:<br>
	> &emsp;&emsp;&emsp;&emsp;arg1 - index of the ST dataset [*options: 1, 2, 3* ] <br>
	> &emsp;&emsp;&emsp;&emsp;arg2 - total number of single-cell reference datasets [*options per removal scenario; <br>
	&emsp;&emsp;&emsp;&emsp;&emsp;rm0= 1:1, rm1= 1:13, rm2= 1:5, rm3= 1:5, rm5= 1:1, rm10= 1:5, rm11= 1:5* ] <br>
	> &emsp;&emsp;&emsp;&emsp;arg3 - removal scenario [*options: rm0, rm1, rm2, rm3, rm5, rm10, rm11* ]
		
		To reproduce the results in the manuscript, the following provides the command line input for no cell type removal (baseline) & one or more cell type removal scenarios,
		Rscript Init_env.R
		(1) ./Execute_R_based_methods.sh 1 1 rm0
		(2) ./Execute_R_based_methods.sh 1 13 rm1
		(3) ./Execute_R_based_methods.sh 1 5 rm2
		(4) ./Execute_R_based_methods.sh 1 5 rm3
		(5) ./Execute_R_based_methods.sh 1 1 rm5
		(6) ./Execute_R_based_methods.sh 1 5 rm10
		(7) ./Execute_R_based_methods.sh 1 5 rm11
		
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; replace `Execute_R_based_methods.sh` by `Execute_python_based_methods.sh` for executing python-based methods. 
	
	*The standalone scripts for each deconvolution method expect the command line arguments as below; notice the different versions of argument* **arg2**.

	> *arg1 - index of the ST dataset <br>
	> arg2 - for R-based methods: index of single-cell reference datasets [varies per removal scenario] <br>
	&emsp;&emsp;&emsp;&emsp;for Python-based methods: total number of single-cell reference datasets <br>
	> arg3 - path to simulated ST datasets directory <br>
	> arg4 - path to single-cell reference datasets directory <br>
	> arg5 - path to save deconvolution results directory* <br>

	Note: While executing cell2location **without GPU support**, `use_gpu` argument should be set to **FALSE** in Cell2Location.py script.

> [!Note]
> In some isolated cases the deconvolution methods indicated below do not generate results, however scripts will still work:<br>
CARD: all 5 instances of the rm11 scenario for all three simulated datasets (ST1, ST2, ST3),<br>
Seurat: for 3 instances of the rm10 scenario for ST1,<br>
Seurat: for 2 instances of the rm11 scenario for ST1.
<br>

### Module 4: Analyse deconvolution results
#### Navigate to "4\_Analysis\_results/Code/" directory 

> [!IMPORTANT]
> Estimating performance metrics needs all the deconvolution results generated in Module 3.
> 
* Estimating performance metrics for all the scenarios of cell type removal.

    - **For JSD calculations**:
   
	&emsp;&emsp;&emsp;&emsp;`Rscript Init_env.R` <br> 
	&emsp;&emsp;&emsp;&emsp;`Rscript Get_JSD_results.R` *arg1 optional\_arg2 optional\_arg3*
	> &emsp;&emsp;&emsp;&emsp;Note: The required command line arguments for `Get_JSD_results.R` script are as below:<br>
	&emsp;&emsp;&emsp;&emsp;arg1 - index of the simulated ST dataset [*options: 1, 2, 3* ] <br>
	&emsp;&emsp;&emsp;&emsp;arg2 - [*optional; default is all scenarios* ] </font> removal scenarios separated by comma (*rm0* is included by default) [*options: rm1, rm2, rm3, rm5, rm10, rm11* ]  <br>
	&emsp;&emsp;&emsp;&emsp;arg3 - [*optional; default is all methods*] name of the deconvolution methods separated by comma <br>
		
		To reproduce the results in the manuscript, the following provides the command line input for calculating JSD metrics for all methods with 1st simulated ST data and all removal of cell type scenarios,
		Rscript Init_env.R
		Rscript Get_JSD_results.R 1

   - **For RMSE calculations**:

	&emsp;&emsp;&emsp;&emsp;`Rscript Get_RMSE_results.R` *arg1 optional\_arg2 optional\_arg3* 
   > &emsp;&emsp;&emsp;&emsp;Note: The required command line arguments for `Get_RMSE_results.R` script are as below:<br>
	&emsp;&emsp;&emsp;&emsp;arg1 - index of the simulated ST dataset [*options: 1, 2, 3* ] <br>
	&emsp;&emsp;&emsp;&emsp;arg2 - [*optional; default is all scenarios* ]  removal scenarios separated by comma (*rm0* is included by default) [*options: rm1, rm2, rm3, rm5, rm10, rm11* ]  <br>
	&emsp;&emsp;&emsp;&emsp;arg3 - [*optional; default is all methods* ] name of the deconvolution methods separated by comma <br>
		
		To reproduce the results in the manuscript, the following provides the command line input for calculating RMSE metrics for all methods with 1st simulated ST data and all removal of cell type scenarios,
		Rscript Get_RMSE_results.R 1

   - **For cell type reassignment calculations:** <br>
   
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;`Rscript Get_celltype_assignment_results.R` *arg1 arg2 optional\_arg3*
	> &emsp;&emsp;&emsp;&emsp;Note: The required command line arguments for `Get_celltype_assignment_results.R` script are as below:<br>
	&emsp;&emsp;&emsp;&emsp;arg1 - index of the simulated ST dataset [*options: 1, 2, 3* ] <br>
	&emsp;&emsp;&emsp;&emsp;arg2 - removal scenarios separated by comma (*rm0* is included by default) [*only supported for rm1,rm2,rm3* ] <br>
	&emsp;&emsp;&emsp;&emsp;arg3 - [*optional; default is all methods* ] name of the deconvolution methods separated by comma <br>
		
		To reproduce the results in the manuscript the following provides the command line input for calculating cell type reassignment for all method with 1st simulated ST data, and removal of one, two & three cell type scenarios,
		Rscript Get_celltype_assignment_results.R 1 "rm1","rm2","rm3"	
> [!Note]
> a list has to be provided as a single argument without spaces in between
<br>

### Module 5: Post analysis results
#### Navigate to "5\_Post\_analysis\_work/Code/" directory 

> [!IMPORTANT]
> This module expects calculation results from Module 4

* Generating figures/plots are included in the manuscript for calculated JSD, RMSE, and cell type reassignment estimates for the specified cell type removal scenarios.
	
   - Generate plots for cell type correlation and cell type reassignment plots for cell type removal scenarios rm1, rm2, and rm3 using commands as below;
   
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;`Rscript Init_env.R` <br> 
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;`Rscript Celltype_correlation.R` <br>
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;`Rscript Fig_celltype_assignments.R` *arg1* <br>
   	> &emsp;&emsp;&emsp;&emsp;Note:`Fig_celltype_assignments.R` expects command line arguments as below:<br>
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;arg1 - index of the simulated ST dataset [*options: 1, 2, 3* ] <br>
	
		To reproduce the results in the manuscript the following provides the command line input,
		Rscript Init_env.R
		Rscript Celltype_correlation.R
		Rscript Fig_celltype_assignments.R 1
		   
   - Generate plots for JSD/RMSE plots for all cell type removal scenarios using commands as below;
   
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;`Rscript Fig_jsd_rmse_plots.R` *arg1* <br>
	> &emsp;&emsp;&emsp;&emsp;Note: `Fig_jsd_rmse_plots.R` expects command line arguments as below; (list of methods and removal scenarios are adapted from JSD/RMSE calculations in Module 4)<br>
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;arg1 - index of the simulated ST dataset [*options: 1, 2, 3* ] <br>
	
		To reproduce the results in the manuscript, the following provides the command line input,
		Rscript Fig_jsd_rmse_plots.R 1

   - Generate a summary funky plot for an overview of the ranking of deconvolution methods across all the removal scenarios using the commands as below;
   
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;`Rscript Funky_plots.R` *arg1* <br>
	> &emsp;&emsp;&emsp;&emsp;Note: `Funky_plots.R` expects command line arguments as below;<br>
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;arg1 - index of the simulated ST dataset [*options: 1, 2, 3* ] <br>
	
		To reproduce the results in the manuscript, the following provides the command line input,
		Rscript Funky_plots.R 1
	
	- For validating the reproduction results and cross-checking the numbers with the original analysis below R script will generate plots that can be compared with those in the original analysis.


	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;`Rscript Validate_numbers.R` *arg1* <br>
	> &emsp;&emsp;&emsp;&emsp;Note: `Validate_numbers.R` expects command line arguments as below;<br>
			&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;arg1 - index of the simulated ST dataset [*options: 1, 2, 3* ] <br>
			
		To reproduce the results in the manuscript, the following provides the command line input,
		Rscript Validate_numbers.R 1
<br>

## Known issues

#### 1. Installation of 'fields' R package on macOS

The analysis is carried out with fields package version 13.3, developed under R version 4.1.2, and RStudio, built for x86_64 architecture.

If you use MacOS with arm\_64 architecture and install RStudio built for x86\_64 architecture, it will use the underlying 'Rosetta2' to run RStudio with x86\_64 architecture, and you will able to install the required 13.3 version of the fields package from CRAN archives.
But if you install RStudio built for arm\_64, it will run with the native silicon chip, and then you will get the error below,


<details>
<summary> <span style="color:blue"> Error message (click to unfold) </span>
</summary>
 
```
Error: Error installing package 'fields':
* installing *source* package ‘fields’ ...
** package ‘fields’ successfully unpacked and MD5 sums checked
** using staged installation
** libs
clang -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/usr/local/include   -fPIC  -Wall -g -O2  -c ExponentialUpperC.c -o ExponentialUpperC.o
clang -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/usr/local/include   -fPIC  -Wall -g -O2  -c RdistEarth.c -o RdistEarth.o
clang -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/usr/local/include   -fPIC  -Wall -g -O2  -c addToDiagC.c -o addToDiagC.o
clang -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/usr/local/include   -fPIC  -Wall -g -O2  -c compactToMatC.c -o compactToMatC.o
clang -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/usr/local/include   -fPIC  -Wall -g -O2  -c expfnC.c -o expfnC.o
/usr/local/bin/gfortran -fno-optimize-sibling-calls  -fPIC  -Wall -g -O2  -c fieldsF77Code.f -o fieldsF77Code.o
fieldsF77Code.f:104:32:

  104 |       double precision A(NMAX,4),V(NMAX,7)
      |                                1
Warning: Array ‘a’ at (1) is larger than limit set by ‘-fmax-stack-var-size=’, moved from stack to static storage. This makes the procedure unsafe when called recursively, or concurrently from multiple threads. Consider increasing the ‘-fmax-stack-var-size=’ limit (or use ‘-frecursive’, which implies unlimited ‘-fmax-stack-var-size’) - or change the code to use an ALLOCATABLE array. If the variable is never accessed concurrently, this warning can be ignored, and the variable could also be declared with the SAVE attribute. [-Wsurprising]
fieldsF77Code.f:108:23:

  108 |       integer idx(NMAX)
      |                       1
Warning: Array ‘idx’ at (1) is larger than limit set by ‘-fmax-stack-var-size=’, moved from stack to static storage. This makes the procedure unsafe when called recursively, or concurrently from multiple threads. Consider increasing the ‘-fmax-stack-var-size=’ limit (or use ‘-frecursive’, which implies unlimited ‘-fmax-stack-var-size’) - or change the code to use an ALLOCATABLE array. If the variable is never accessed concurrently, this warning can be ignored, and the variable could also be declared with the SAVE attribute. [-Wsurprising]
fieldsF77Code.f:107:23:

  107 |       integer imx(NMAX)
      |                       1
Warning: Array ‘imx’ at (1) is larger than limit set by ‘-fmax-stack-var-size=’, moved from stack to static storage. This makes the procedure unsafe when called recursively, or concurrently from multiple threads. Consider increasing the ‘-fmax-stack-var-size=’ limit (or use ‘-frecursive’, which implies unlimited ‘-fmax-stack-var-size’) - or change the code to use an ALLOCATABLE array. If the variable is never accessed concurrently, this warning can be ignored, and the variable could also be declared with the SAVE attribute. [-Wsurprising]
fieldsF77Code.f:106:59:

  106 |       double precision ux(NMAX),uy(NMAX), uw(NMAX),ud(NMAX),utr
      |                                                           1
Warning: Array ‘ud’ at (1) is larger than limit set by ‘-fmax-stack-var-size=’, moved from stack to static storage. This makes the procedure unsafe when called recursively, or concurrently from multiple threads. Consider increasing the ‘-fmax-stack-var-size=’ limit (or use ‘-frecursive’, which implies unlimited ‘-fmax-stack-var-size’) - or change the code to use an ALLOCATABLE array. If the variable is never accessed concurrently, this warning can be ignored, and the variable could also be declared with the SAVE attribute. [-Wsurprising]
fieldsF77Code.f:106:50:

  106 |       double precision ux(NMAX),uy(NMAX), uw(NMAX),ud(NMAX),utr
      |                                                  1
Warning: Array ‘uw’ at (1) is larger than limit set by ‘-fmax-stack-var-size=’, moved from stack to static storage. This makes the procedure unsafe when called recursively, or concurrently from multiple threads. Consider increasing the ‘-fmax-stack-var-size=’ limit (or use ‘-frecursive’, which implies unlimited ‘-fmax-stack-var-size’) - or change the code to use an ALLOCATABLE array. If the variable is never accessed concurrently, this warning can be ignored, and the variable could also be declared with the SAVE attribute. [-Wsurprising]
fieldsF77Code.f:106:31:

  106 |       double precision ux(NMAX),uy(NMAX), uw(NMAX),ud(NMAX),utr
      |                               1
Warning: Array ‘ux’ at (1) is larger than limit set by ‘-fmax-stack-var-size=’, moved from stack to static storage. This makes the procedure unsafe when called recursively, or concurrently from multiple threads. Consider increasing the ‘-fmax-stack-var-size=’ limit (or use ‘-frecursive’, which implies unlimited ‘-fmax-stack-var-size’) - or change the code to use an ALLOCATABLE array. If the variable is never accessed concurrently, this warning can be ignored, and the variable could also be declared with the SAVE attribute. [-Wsurprising]
fieldsF77Code.f:106:40:

  106 |       double precision ux(NMAX),uy(NMAX), uw(NMAX),ud(NMAX),utr
      |                                        1
Warning: Array ‘uy’ at (1) is larger than limit set by ‘-fmax-stack-var-size=’, moved from stack to static storage. This makes the procedure unsafe when called recursively, or concurrently from multiple threads. Consider increasing the ‘-fmax-stack-var-size=’ limit (or use ‘-frecursive’, which implies unlimited ‘-fmax-stack-var-size’) - or change the code to use an ALLOCATABLE array. If the variable is never accessed concurrently, this warning can be ignored, and the variable could also be declared with the SAVE attribute. [-Wsurprising]
fieldsF77Code.f:104:42:

  104 |       double precision A(NMAX,4),V(NMAX,7)
      |                                          1
Warning: Array ‘v’ at (1) is larger than limit set by ‘-fmax-stack-var-size=’, moved from stack to static storage. This makes the procedure unsafe when called recursively, or concurrently from multiple threads. Consider increasing the ‘-fmax-stack-var-size=’ limit (or use ‘-frecursive’, which implies unlimited ‘-fmax-stack-var-size’) - or change the code to use an ALLOCATABLE array. If the variable is never accessed concurrently, this warning can be ignored, and the variable could also be declared with the SAVE attribute. [-Wsurprising]
fieldsF77Code.f:379:43:

  379 |       double precision work(nobs),diag(mxM),dumm1(1),dumm2(1)
      |                                           1
Warning: Array ‘diag’ at (1) is larger than limit set by ‘-fmax-stack-var-size=’, moved from stack to static storage. This makes the procedure unsafe when called recursively, or concurrently from multiple threads. Consider increasing the ‘-fmax-stack-var-size=’ limit (or use ‘-frecursive’, which implies unlimited ‘-fmax-stack-var-size’) - or change the code to use an ALLOCATABLE array. If the variable is never accessed concurrently, this warning can be ignored, and the variable could also be declared with the SAVE attribute. [-Wsurprising]
clang -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/usr/local/include   -fPIC  -Wall -g -O2  -c init.c -o init.o
clang -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/usr/local/include   -fPIC  -Wall -g -O2  -c multebC.c -o multebC.o
clang -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/usr/local/include   -fPIC  -Wall -g -O2  -c rdistC.c -o rdistC.o
clang -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -L/Library/Frameworks/R.framework/Resources/lib -L/usr/local/lib -o fields.so ExponentialUpperC.o RdistEarth.o addToDiagC.o compactToMatC.o expfnC.o fieldsF77Code.o init.o multebC.o rdistC.o -L/usr/local/gfortran/lib/gcc/x86_64-apple-darwin18/8.2.0 -L/usr/local/gfortran/lib -lgfortran -lquadmath -lm -F/Library/Frameworks/R.framework/.. -framework R -Wl,-framework -Wl,CoreFoundation
ld: warning: -single_module is obsolete
ld: warning: -multiply_defined is obsolete
ld: warning: search path '/usr/local/gfortran/lib/gcc/x86_64-apple-darwin18/8.2.0' not found
ld: warning: ignoring file '/private/var/folders/16/2qzgz8cd2bl65zhgv1x6sy0c0000gn/T/RtmpcY8v1T/R.INSTALL3349486d5ed2/fields/src/fieldsF77Code.o': found architecture 'arm64', required architecture 'x86_64'
ld: warning: ignoring file '/usr/local/lib/libgfortran.5.dylib': found architecture 'arm64', required architecture 'x86_64'
ld: warning: ignoring file '/usr/local/gfortran/lib/libquadmath.0.dylib': found architecture 'arm64', required architecture 'x86_64'
installing to /Users/utkarsh/surfdrive/UtkarshMahamune/20210701_RobustnessEvaluation/Processing/renv/staging/1/00LOCK-fields/00new/fields/libs
** R
** data
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
** building package indices
** testing if installed package can be loaded from temporary location
Error: package or namespace load failed for ‘fields’ in dyn.load(file, DLLpath = DLLpath, ...):
 unable to load shared object '/Users/utkarsh/surfdrive/UtkarshMahamune/20210701_RobustnessEvaluation/Processing/renv/staging/1/00LOCK-fields/00new/fields/libs/fields.so':
  dlopen(/Users/utkarsh/surfdrive/UtkarshMahamune/20210701_RobustnessEvaluation/Processing/renv/staging/1/00LOCK-fields/00new/fields/libs/fields.so, 0x0006): symbol not found in flat namespace '_css_'
Error: loading failed
Execution halted
ERROR: loading failed
* removing ‘/Users/utkarsh/surfdrive/UtkarshMahamune/20210701_RobustnessEvaluation/Processing/renv/staging/1/fields’
install of package 'fields' failed [error code 1]
Traceback (most recent calls last):
14: renv::init()
13: restore(project = project, library = libpaths, repos = repos, 
        prompt = FALSE)
12: renv_restore_run_actions(project, diff, current, lockfile, rebuild)
11: renv_install_impl(records)
10: renv_install_staged(records)
 9: renv_install_default(records)
 8: handler(package, renv_install_package(record))
 7: renv_install_package(record)
 6: withCallingHandlers(renv_install_package_impl(record), error = function(e) writef("FAILED"))
 5: renv_install_package_impl(record)
 4: r_cmd_install(package, path)
 3: r_exec_error(package, output, "install", status)
 2: abort(all)
 1: stop(fallback)

```
 
</details>

to resolve this error, you can follow either option from below, 

- You can install the latest version 15.2 of fields R package ([R-binaries](https://cran.r-project.org/web/packages/fields/index.html) available for arm\_64); this will affect the CARD package since it needs to be updated to the latest version as well.

- Install R and RStudio built for x86\_64 architecture and run the analysis using it. If you still get the same error, the compiler uses the default arm\_64 architecture to install R packages that need compilation.

<br>

#### 2. RCTD parallel execution

If you come across the below error while executing the RCTD deconvolution method.

```
Error in checkForRemoteErrors(lapply(cl, recvResult)) :
  4 nodes produced errors; first error: object '.doSnowGlobals' not found
Calls: run.RCTD ... %dopar% -> <Anonymous> -> clusterCall -> checkForRemoteErrors
Execution halted
EEError in unserialize(node$con) : error reading from connection
Calls: <Anonymous> ... doTryCatch -> recvData -> recvData.SOCKnode -> unserialize
rror in unserialize(node$con) : error reading from connection
Calls: <Anonymous> ... doTryCatch -> recvData -> recvData.SOCKnode -> unserialize
EError in unserialize(node$con) : error reading from connection
Calls: <Anonymous> ... doTryCatch -> recvData -> recvData.SOCKnode -> unserialize
Error in unserialize(node$con) : error reading from connection
Calls: <Anonymous> ... doTryCatch -> recvData -> recvData.SOCKnode -> unserialize
Execution halted
Execution halted
Execution halted
Execution halted
```

This is a known issue when RCTD runs multiple jobs in parallel from within the deconvolution function. This has been reported to the developers earlier and can be found in the issues on the GitHub repository ([link](https://github.com/dmcable/spacexr/issues/141)).

Please follow the solution available in the reported issue; if the problem persists, please raise an issue on GitHub.

<br>

#### 3. Creating conda virtual environment with environment.yml


The `environment.yml` is generated on the linux system and thus you can get error like below.

```
PackagesNotFoundError: The following packages are not available from current channels:

  - libstdcxx-ng=11.2.0*
  - libgomp=11.2.0*
  - libgcc-ng=11.2.0*
  - ld_impl_linux-64=2.40*
  - _openmp_mutex=5.1*
  - _libgcc_mutex=0.1*
```
<br>

#### 4. python\_config\_impl(python) error

The R script uses installed Python using the RETICULATE package. You get this error if the package cannot find the proper Python installation.

The issue can be resolved by adding the line below at the beginning of the R script that is giving the error to set up the path for RETICULATE_python.

`Sys.setenv(RETICULATE_PYTHON="path_of_the_desired_python_installation")`

<br>

## Appendix A
### Installation steps for Windows OS

1. Install R 4.1.2 [download here](https://cran.r-project.org/bin/windows/base/old/4.1.2/R-4.1.2-win.exe)
2. Install RTools 4.0  [download here](https://cran.r-project.org/bin/windows/Rtools/rtools40.html). <br>
Add RTools path to *'PATH'* environment variable (make sure no white spaces are present in the path)
3. Install miniconda [download here](https://repo.anaconda.com/miniconda/) <br>
Look for *'Miniconda3-py39_24.11.1-0-Windows-x86_64.exe'*
4. Set up conda environment with Python version 3.9.7 using the command: <t>`conda env create -f environment.yml`
5. Set up renv setup details
	- Execute `Renv_setup.R` script from where it resides currently to initialise the renv infrastructure for the R project (ignore the warning messages) using `Rscript Renv_setup.R arg1` command. <br>
	Expected command line argument (*arg1*) is **(classic) GITHUB_PAT token** ([more details](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens)).<br>
	*The required files for renv setup (renv.lock, .Rprofile, renv/activate.R, renv/settings.json) should already be in your R project directory; If not, ensure it resides in the same directory as the Renv_setup.R file.*
