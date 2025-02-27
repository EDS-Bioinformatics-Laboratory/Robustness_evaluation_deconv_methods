# <font color="darkgreen"> Systematic evaluation of robustness to cell type mismatch of deconvolution methods for spatial transcriptomics data </font>

&emsp;&emsp;&emsp;&emsp; Spatial transcriptomics approaches based on sequencing (barcode-based, e.g., 10x Visium) preserve spatial information but with limited cellular resolution. On the other hand, single-cell RNA-sequencing (scRNA-seq) techniques provide single-cell resolution but lose spatial resolution because of the tissue dissociation step during the scRNA-seq experimental procedure. With these complementary strengths in mind, computational methods have been developed to combine scRNA-seq and spatial transcriptomics data. These approaches use deconvolution to identify cell types and their respective proportions at each location in spatial transcriptomics data with the aid of a scRNA-seq reference dataset. Some suggest that deconvolution approaches are sensitive to the absence of cell type(s) in the single-cell reference dataset, a problem referred to as *cell type mismatch*.

Here, we systematically evaluated the robustness of deconvolution methods to cell type mismatch tailored for spatial transcriptomics data.
<br>
<br>

## <font color = "green">Content of the Processing directory</font>

* **0_SoftwareEnvironment**:
This directory enlists the software environment specifications used for various programming languages and/or platforms during the project.

* **1\_Generate\_sc\_ref\_data**: The directory comprises scripts, results and settings to generate the integrated scRNA-seq dataset from 2 distinct and complementary scRNA-seq datasets. The final result **`sc.ref.data.rds`** is used as a basis to generate various reference dataset versions based on cell type removal scenarios and to generate simulated spatial transcriptomics datasets.

* **2\_Simulating\_ST\_data**: Comprised of scripts, results and settings to generate the sequencing-based simulated spatial transcriptomics datasets varying in number of cells & cell types present per spatial location (spot). We created three simulated ST datasets using the algorithm developed for the analysis. <br>
&emsp;&emsp;&emsp; a. **<font color="brown">ST1</font>**: 4-8 cell types and 10-15 cells per spot; <br>
&emsp;&emsp;&emsp; b. **<font color="brown">ST2</font>**: 1-5 cell types and 10-15 cells per spot; <br>
&emsp;&emsp;&emsp; c. **<font color="brown">ST3</font>**: 1-5 cell types and 3-7 cells per spot.


* **3\_ST\_methods**: Comprised of standalone R/Python scripts, one for each deconvolution method and shell/batch scripts to execute the R or Python scripts in parallel for multiple instances in a removal scenario, <font color="blue">provided the required computational power is available</font>. Six of the eight methods are *R-based* (CARD, RCTD, SCDC, MuSiC, Seurat, SPOTlight), while two are *Python-based* (cell2location, Stereoscope). The Python-based methods <font color="blue"> expect GPU support </font>.<br>
&emsp;&emsp;&emsp; The resulting deconvolution results for each removal scenario are available in the respective directory (for instance, the _**rm1**_ directory refers to the removal scenario for removing one cell type from scRNA-seq reference data). See below the overview of removal scenarios and total reference datasets in each scenario. <br>
&emsp;&emsp;&emsp; 1. **<font color="brown">rm0</font>** - removal of _**no**_ cell types from reference data | one reference dataset <br>
&emsp;&emsp;&emsp; 2. **<font color="brown">rm1</font>** - removal of _**one**_ cell type from reference data | 13 reference datasets <br>
&emsp;&emsp;&emsp; 3. **<font color="brown">rm2</font>** - removal of _**two**_ cell types from reference data | 5 reference datasets <br>
&emsp;&emsp;&emsp; 4. **<font color="brown">rm3</font>** - removal of _**three**_ cell types from reference data | 5 reference datasets <br>
&emsp;&emsp;&emsp; 5. **<font color="brown">rm5</font>** - removal of _**five**_ cell types from reference data | 1 reference dataset <br>
&emsp;&emsp;&emsp; 6. **<font color="brown">rm10</font>** - removal of _**ten**_ cell types from reference data | 5 reference datasets <br>
&emsp;&emsp;&emsp; 7. **<font color="brown">rm11</font>** - removal of _**eleven**_ cell types from reference data | 5 reference datasets


* **4\_Analysis\_results**: Comprised of scripts to calculate similarity metrics like JSD and RMSE to understand the performance of deconvolution methods for various cell type removal scenarios compared to baseline with no cell type missing. Cell type reassignment metrics calculate the assignment of missing cell type proportions from the reference dataset.

* **5\_Post\_analysis\_work**: The directory comprises scripts to generate the final resulting plots in the project; few plots are included in the manuscript's main text, while others were included in the supplementary information.

* **Data**: Comprised of pre-processed data downloaded from a public data-sharing platform. One of the datasets is procured from the group of Lisa van Baarsen and is available on GEO as well with GEO accession ID - ####. <font color="brown">(not available in github repository)</font>

* **renv**: The project uses the renv functionality and the directory comprises the essential files and directories for the renv setup.

* **Renv_setup.R**: This script sets up the R environment for the project work. Details about executing are available in the file as a header and in the comments. Also included under the 'How to reproduce the results section further dowm in this document.<br>
**Please read the instructions in the file carefully before proceeding if you are using RStudio.**

* **renv.lock**: R environment lock (metadata) file comprising package details. This file will be used by the <font color="brown">*Renv_setup.R*</font> script to download and install the correct package and its version to recreate the computing environment.

* **.Rprofile**: Essential R profile script to set up the correct R profile while using the renv infrastructure.

* **environment.yml**: A metadata file comprising of python and/or packages installed in the virtual conda environment.

<br>

## <font color = "green"> How can you reproduce the results in the manuscript? </font>


<font color="red" size = 3> **Note**: The analysis has been carried out using R version <u>**4.1.2**</u> [tested on macos/ubuntu/windows platform], and Python version <u>**3.9.7**</u> [tested on ubuntu platform]; please make sure you use the same version while trying to reproduce the results</font>

* The software environment for R and Python used during the original analysis can be found under the `/0_SoftwareEnvironment` directory.


* For reproducibility, the project utilises conda and renv functionality. <br>

	<font color="darkblue" size=4>§ Steps to set up conda functionalities for your analysis </font> <br>
	
	**I**. 	<font color="brown">[optional]</font> Read [<font color="blue">getting started with conda</font>](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html) carefully to understand the functionalities provided by the conda environment.
	
	**II.**	If you do not already have conda installed on your machine, please install it by following [<font color="blue"> steps to install conda</font>](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).
	
	**III**. To create a conda virtual environment similar to the analysis, use the <font color="brown">*environment.yml*</font> file available in the root directory or under `0_SoftwareEnvironment/Python` directory. The command for creating a conda virtual environment from the command line terminal/prompt is as below;
	
	&emsp;&emsp;&emsp;&emsp;<font color="blue" size=4>`conda env create -f environment.yml`</font> 

	> **<font color="red">Note:</font>** Depending on your platform(OS), installation of a few Python packages/libraries might fail; please check the logs/status of the environment creation carefully.<br>

	<font color="darkblue" size=4>§ Steps to set up renv functionalities for your analysis </font> <br>
	
	**I**. 	<font color="brown">[optional]</font> Read [<font color="blue">getting started with renv</font>](https://rstudio.github.io/renv/articles/renv.html) carefully to understand the functionalities provided by the renv R package.
		
	**II**. Executing `Renv_Setup.R` from where it resides currently to initialise the renv infrastructure for the R project (ignore the warning messages):
	- The `Renv_Setup.R` script expects a **GITHUB_PAT token (classic)** as command line argument (arg1). More details on the github token can be found [here](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens).<br>

		&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue" size=4>`Rscript Renv_setup.R arg1` </font> 
		
	&emsp;&emsp;&emsp;&emsp;
*The required files for renv setup <font color="brown">(renv.lock, .Rprofile, renv/activate.R, renv/setting.json)</font> should already be in your R project directory; If not, ensure it resides in the same directory as the Renv_setup.R file.*
	
	> **<font color="red">Note:</font>** A collaborator can still experience minor discrepancies in the results while using the renv functionality due to the platform(OS) dependencies.
	
<br>
<br>

**<font color="red" size=4>Disclaimer:</font>** 
<font color="red">The following instructions are for reproducing the results using the command-line terminal/prompt, not RStudio or Jupyter Notebook. </font>
<br>
<br>

#### <font color="darkblue">Module-1: Navigate to "1\_Generate\_sc\_ref_data/Code/" directory </font>

* Generating a single-cell reference dataset from multiple single-cell datasets, UMAP representations included in the supplementary information of the manuscript using the commands as below; 
	
	&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Init_env.R`</font> <br>
	&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript SC_ref_data.R` </font> 
	
* Generating single-cell reference datasets for all removal scenarios based on cell type removal using the commands as below;
	
	&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Init_env.R`</font> <br>
	&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript SC_ref_data_scenarios.R`</font>
	> <font color="red">Note:</font>`SC_ref_data_scenarios.R` script requires the single-cell reference dataset generated in the previous step.
	
<br>

#### <font color="darkblue">Module-2: Navigate to "2\_Simulating\_ST\_data/Code/" directory </font>

* Generating simulated spatial transcriptomics datasets using single-cell reference data, the three datasets vary by the number of cells and cell types present per spatial location using the commands as below;

	&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Init_env.R`</font> <br>
	&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Generate_ST_data.R` *arg1 arg2 arg3 arg4 arg5*</font>
	
	> &emsp;&emsp;&emsp;<font color="red">Note:</font> `Generate_ST_data.R` script expects five command line arguments in the order mentioned below;<br>
	>	&emsp;&emsp;&emsp;&emsp;arg1 - min number of cell types per spot <br>
	>	&emsp;&emsp;&emsp;&emsp;arg2 - max number of cell types per spot <br>
	>	&emsp;&emsp;&emsp;&emsp;arg3 - min number of cells per spot <br>
	>	&emsp;&emsp;&emsp;&emsp;arg4 - max number of cells per spot <br>
	>	&emsp;&emsp;&emsp;&emsp;arg5 - index of simulated ST dataset [*options: 1, 2, 3* ] <br>
	
		command line input to generate 1st simulated ST dataset as below,
		Rscript Generate_ST_data.R 4 8 10 15 1

<br>

#### <font color="darkblue">Module-3: Navigate to "3\_ST\_methods/Code/" directory </font>

* The shell/batch scripts execute all deconvolution methods to predict cell type proportions simultaneously for all removal scenarios and multiple reference datasets for each scenario. <br>

	Shell/batch scripts to execute R/Python-based deconvolution methods: <br>
		<font color="blue">`Rscript Init_env.R`</font> <br>
	<font color="blue">`Execute_R_based_methods.sh` *arg1 arg2 arg3*</font>, <br>
	<font color="blue">`Execute_python_based_methods.sh` *arg1 arg2 arg3*</font>, <br>
	<br>
		<font color="blue">`Execute_R_based_methods.bat` *arg1 arg2 arg3*</font>, <br>
	<font color="blue">`Execute_python_based_methods.bat` *arg1 arg2 arg3*</font>
	
	Please make sure to execute only one script at a time and wait for the results. For every removal scenario, each single-cell reference will have one result (for instance, rm1 will have 13 results for each method) for every method. Verify the number of results generated moving to the next scenario. If a combination of SC & ST dataset is failed for a method, please run the standalone script for that particular pair. <br>

	> <font color="red">Note:</font> The script to execute R-based/Python-based methods executes each method in parallel for a provided number of reference datasets and removal scenarios; the required command line arguments are as below:<br>
	> &emsp;&emsp;&emsp;&emsp;arg1 - index of the ST dataset [*options: 1, 2, 3* ] <br>
	> &emsp;&emsp;&emsp;&emsp;arg2 - total number of single-cell reference datasets [*options per removal scenario; <br>
	&emsp;&emsp;&emsp;&emsp;&emsp;rm0= 1:1, rm1= 1:13, rm2= 1:5, rm3= 1:5, rm5= 1:1, rm10= 1:5, rm11= 1:5* ] <br>
	> &emsp;&emsp;&emsp;&emsp;arg3 - removal scenario [*options: rm0, rm1, rm2, rm3, rm5, rm10, rm11* ]
		
		command line input for no cell type removal (baseline) & one or more cell type removal scenarios would be then as follows,
		(1) ./Execute_R_based_methods.sh 1 1 rm0
		(2) ./Execute_R_based_methods.sh 1 13 rm1
		(3) ./Execute_R_based_methods.sh 1 5 rm2
		(4) ./Execute_R_based_methods.sh 1 5 rm3
		(5) ./Execute_R_based_methods.sh 1 1 rm5
		(6) ./Execute_R_based_methods.sh 1 5 rm10
		(7) ./Execute_R_based_methods.sh 1 5 rm11
		

	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; replace `Execute_R_based_methods.sh` by `Execute_python_based_methods.sh` for executing python-based methods. 
	
	*The standalone scripts for each deconvolution method expect the command line arguments as below; notice the different versions of argument <font color="red">b </font>*.

	> *arg1 - index of the ST dataset <br>
	> arg2 - for R-based methods: index of single-cell reference datasets [varies per removal scenario] <br>
	&emsp;&emsp;&emsp;&emsp;for Python-based methods: total number of single-cell reference datasets <br>
	> arg3 - path to simulated ST datasets directory <br>
	> arg4 - path to single-cell reference datasets directory <br>
	> arg5 - path to save deconvolution results directory* <br>
	

	<font color="red"> **Note**: While executing cell2location **without GPU support**, `use_gpu` argument should be set to **FALSE** in Cell2Location.py script.</font>

<br>

#### <font color="darkblue">Module-4: Navigate to "4\_Analysis\_results/Code/" directory </font>

* Estimating performance metrics for all the removal of cell types scenarios.

	> <font color="red">Note:</font> Estimating performance metrics needs all the deconvolution results generated in Module-3.

    - **For JSD calculations**:
   
	&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Init_env.R`</font> <br> 
	&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Get_JSD_results.R` *arg1 optional\_arg2 optional\_arg3*</font>
	> &emsp;&emsp;&emsp;&emsp;<font color="red">Note: </font> The required command line arguments for <font color="blue">`Get_JSD_results.R`</font>script are as below:<br>
	&emsp;&emsp;&emsp;&emsp;arg1 - index of the simulated ST dataset [*options: 1, 2, 3* ] <br>
	&emsp;&emsp;&emsp;&emsp;arg2 - [*optional; default is all scenarios* ] </font> removal scenarios separated by comma (*rm0* is included by default) [*options: rm1, rm2, rm3, rm5, rm10, rm11* ]  <br>
	&emsp;&emsp;&emsp;&emsp;arg3 - [*optional; default is all methods*] name of the deconvolution methods separated by comma <br>
		
		command line input for calculating JSD metrics for all methods with 1st simulated ST data and all removal of cell type scenarios would be,
		Rscript Get_JSD_results.R 1

   - **For RMSE calculations**:

   &emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Init_env.R` </font> <br> 
	&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Get_RMSE_results.R` *arg1 optional\_arg2 optional\_arg3* </font>
   > &emsp;&emsp;&emsp;&emsp;<font color="red">Note: </font> The required command line arguments for <font color="blue">`Get_RMSE_results.R`</font>script are as below:<br>
	&emsp;&emsp;&emsp;&emsp;arg1 - index of the simulated ST dataset [*options: 1, 2, 3* ] <br>
	&emsp;&emsp;&emsp;&emsp;arg2 - [*optional; default is all scenarios* ]  removal scenarios separated by comma (*rm0* is included by default) [*options: rm1, rm2, rm3, rm5, rm10, rm11* ]  <br>
	&emsp;&emsp;&emsp;&emsp;arg3 - [*optional; default is all methods* ] name of the deconvolution methods separated by comma <br>
		
		command line input for calculating RMSE metrics for all methods with 1st simulated ST data and all removal of cell type scenarios would be,
		Rscript Get_RMSE_results.R 1

   - **For cell type reassignment calculations:** <br>
   
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Init_env.R`</font> <br> 
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Get_celltype_assignment_results.R` *arg1 arg2 optional\_arg3*</font>
	> &emsp;&emsp;&emsp;&emsp;<font color="red">Note: </font> The required command line arguments for <font color="blue">`Get_celltype_assignment_results.R`</font> script are as below:<br>
	&emsp;&emsp;&emsp;&emsp;arg1 - index of the simulated ST dataset [*options: 1, 2, 3* ] <br>
	&emsp;&emsp;&emsp;&emsp;arg2 - removal scenarios separated by comma (*rm0* is included by default) [*options: rm1, rm2, rm3* ] <br>
	&emsp;&emsp;&emsp;&emsp;arg3 - [*optional; default is all methods* ] name of the deconvolution methods separated by comma <br>
		
		command line input for calculating cell type reassignment for method CARD and RCTD with 1st simulated ST data, and removal of one, two & three cell type scenarios would be,
		Rscript Get_celltype_assignment_results.R 1 "rm1","rm2","rm3" "CARD","RCTD"
	
	Note that a list is provided as a single argument without spaces in between

<br>

#### <font color="darkblue">Module-5: Navigate to "5\_Post\_analysis\_work/Code/" directory </font>

* Generating figures/plots included in the manuscript for calculated JSD, RMSE and cell type reassignment estimates for the specified cell type removal scenarios.

	> <font color="red">Note:</font> this module expects calculations results from Module-4
	
   - Generate plots for cell type correlation and cell type reassignment plots for cell type removal scenarios rm1, rm2, and rm3 using commands as below;
   
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Init_env.R`</font> <br> 
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Celltype_correlation.R `</font> <br>
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Fig_celltype_assignments.R` *arg1*</font> <br>
   	> &emsp;&emsp;&emsp;&emsp;<font color="red">Note:</font> <font color="blue">`Fig_celltype_assignments.R`</font> expects command line arguments as below:<br>
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;arg1 - index of the simulated ST dataset [*options: 1, 2, 3* ] <br>
   
   - Generate plots for JSD/RMSE plots for all cell type removal scenarios using commands as below;
   
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Init_env.R`</font> <br> 
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Fig_jsd_rmse_plots.R` *arg1*</font> <br>
	> &emsp;&emsp;&emsp;&emsp;<font color="red">Note:</font> <font color="blue">`Fig_jsd_rmse_plots.R`</font> expects command line arguments as below; (<font color="brown" size=2>list of methods and removal scenarios are adapted from JSD/RMSE calculations in Module-4 </font>)<br>
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;arg1 - index of the simulated ST dataset [*options: 1, 2, 3* ] <br>

   - Generate summary funky plot for an overview of ranking of deconvolution methods across all the removal scenarios using commands as below;
   
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Init_env.R`</font> <br> 
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Funky_plots.R` *arg1*</font> <br>
	> &emsp;&emsp;&emsp;&emsp;<font color="red">Note:</font> <font color="blue">`Funky_plots.R`</font> expects command line arguments as below;<br>
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;arg1 - index of the simulated ST dataset [*options: 1, 2, 3* ] <br>
	
	- For validating the reproduction results and cross-checking the numbers with the original analysis below R script will generate plots that can be compared with those in the original analysis.

	   
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Init_env.R`</font> <br> 
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Validate_numbers.R` *arg1*</font> <br>
	> &emsp;&emsp;&emsp;&emsp;<font color="red">Note:</font> <font color="blue">`Validate_numbers.R `</font> expects command line arguments as below;<br>
			&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;arg1 - index of the simulated ST dataset [*options: 1, 2, 3* ] <br>
			

<br>
<br>

## <font color = "red"> Known issues </font>

#### 1. Installation of 'fields' R package on macOS

The analysis is carried out with fields package version 13.3 developed under R version 4.1.2 and RStudio built for x86\_64 architecture.

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

#### 2. Compilation failed error on Windows OS

R and RStudio on Windows OS require rtools (a toolchain for building R and R packages). The correct version of rtools can be downloaded from [here](https://cran.r-project.org/bin/windows/Rtools/).

> <font color="red">Note:</font> Ensure you do not have white spaces in your paths while installing the rtools.
	
<br>

#### 3. Miniconda installation prompt

		No non-system installation of Python could be found.
		Would you like to download and install Miniconda?
		Miniconda is an open-source environment management system for Python.
		See https://docs.conda.io/en/latest/miniconda.html for more details.
		 
		Would you like to install Miniconda? [Y/n]:

If you see the above prompt in your execution attempt, press <font color="red">**n**</font> and stop the execution.

Please check your Python installation before proceeding. You can check this using `Sys.which("python")` or `Sys.which("python3")` in the RStudio console.

- If you see `" "`, please install the python version 3.9.7 (*or higher*) and try again.
- If you see `"/path_to_python_installation"`, please contact the developer or create an issue on the GitHub repository with the error message, R session info, and Python installation details (path, version, etc).

<br>

#### 4. RCTD parallel execution

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

#### 5. Creating conda virtual environment with environment.yml


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


Please comment out the lines in the `environment.yml` file with the package names that gave the error and rerun the command.
