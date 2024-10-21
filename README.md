<<<<<<< HEAD
# <font color="darkgreen"> Systematic evaluation of robustness to cell type mismatch of deconvolution methods for spatial transcriptomics data </font>

&emsp;&emsp;&emsp;&emsp; Spatial transcriptomics approaches based on sequencing (barcode-based, e.g., 10x Visium) preserve spatial information but with limited cellular resolution. On the other hand, single-cell RNA-sequencing (scRNA-seq) techniques provide single-cell resolution but lose spatial resolution because of the tissue dissociation step during the scRNA-seq experimental procedure. With these complementary strengths in mind, computational methods have been developed to combine scRNA-seq and spatial transcriptomics data. These approaches use deconvolution to identify cell types and their respective proportions at each location in spatial transcriptomics data with the aid of a scRNA-seq reference dataset. Some suggest that deconvolution approaches are sensitive to the absence of cell type(s) in the single-cell reference dataset, a problem referred to as *cell type mismatch*.
=======
### Systematic evaluation of robustness to cell type mismatch of cell type deconvolution methods for spatial transcriptomics data 
>>>>>>> parent of c830ac4 (Updated scripts for better reproducibility)

Spatial transcriptomics approaches based on sequencing (barcode-based, e.g., 10x Visium) preserve spatial information but with limited cellular resolution. On the other hand, single-cell RNA-sequencing (scRNA-seq) techniques provide single-cell resolution but lose spatial resolution because of the tissue dissociation step during the scRNA-seq protocol. With these complementary strengths in mind, computational methods have been developed to combine scRNA-seq and spatial transcriptomics data. These approaches use deconvolution to identify cell types and their respective proportions at each location in spatial transcriptomics data with the aid of a scRNA-seq reference dataset. Some suggest that deconvolution approaches are sensitive to the absence of cell type(s) in the single-cell reference dataset, a problem referred to as *cell type mismatch*.
Here, we systematically evaluated the robustness of deconvolution methods to cell type mismatch tailored for spatial transcriptomics data.

#### Content of the Processing directory

**0_SoftwareEnvironment**:
This directory enlists the software environment specifications used during the project.

**Data**: Comprised of pre-processed data either downloaded from a public data-sharing platform or procured from the group of Lisa van Baarsen.

**1\_Generate\_sc\_ref\_data**: The directory comprises scripts, results and settings to generate the integrated scRNA-seq dataset from 2 distinct and complementary scRNA-seq datasets. The final result `sc.ref.data.rds` is used in the further downstream analysis as a reference dataset for generating simulated spatial transcriptomics datasets and as an input to deconvolution methods.

<<<<<<< HEAD
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; a. **<font color="brown">ST1</font>**: 4-8 cell types and 10-15 cells per spot; <br>
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; b. **<font color="brown">ST2</font>**: 1-5 cell types and 10-15 cells per spot; <br>
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; c. **<font color="brown">ST3</font>**: 1-5 cell types and 3-7 cells per spot.
=======
**2\_Simulating\_ST\_data**: Comprised of scripts, results and settings to generate the sequencing-based simulated spatial transcriptomics datasets varying in number of cells & cell types present per spatial location. We created 4 simulated ST datasets using the algorithm developed for the analysis.
>>>>>>> parent of c830ac4 (Updated scripts for better reproducibility)

**3\_ST\_methods**: Comprised of standalone R/python scripts, one for each deconvolution method and shell scripts to execute the R or python scripts in parallel, provided the required computational power is available. The resulting deconvolution results for each removal scenario are available in respected directory (for instance rm1 directory referes to removal of 1 cell type from scRNA-seq reference data).

<<<<<<< HEAD
* **3\_ST\_methods**: Comprised of standalone R/Python scripts, one for each deconvolution method and shell/batch scripts to execute the R or Python scripts in parallel for multiple instances in a removal scenario, <font color="blue">provided the required computational power is available</font>. Six of the eight methods are R-based, while two are Python-based. The Python-based methods <font color="blue"> expect GPU support </font>.<br>
&emsp;&emsp;&emsp;&emsp;The resulting deconvolution results for each removal scenario are available in the respected directory (for instance, the _**rm1**_ directory refers to the removal scenario for removing one cell type from scRNA-seq reference data). See below the overview of removal scenarios and total reference datasets in each scenario. <br>
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; 1. **<font color="brown">rm0</font>** - removal of _**no**_ cell types from reference data &emsp;| one reference dataset <br>
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; 2. **<font color="brown">rm1</font>** - removal of _**one**_ cell type from reference data &emsp;| 13 reference datasets <br>
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; 3. **<font color="brown">rm2</font>** - removal of _**two**_ cell types from reference data &emsp;| 5 reference datasets <br>
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; 4. **<font color="brown">rm3</font>** - removal of _**three**_ cell types from reference data &emsp;| 5 reference datasets <br>
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; 5. **<font color="brown">rm5</font>** - removal of _**five**_ cell types from reference data &emsp;| 1 reference dataset <br>
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; 6. **<font color="brown">rm10</font>** - removal of _**ten**_ cell types from reference data &emsp;| 5 reference datasets <br>
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; 7. **<font color="brown">rm11</font>** - removal of _**eleven**_ cell types from reference data &emsp;| 5 reference datasets
=======
**4\_Analysis\_results**: Comprised of scripts to calculate similarity metrics like JSD and RMSE to understand the performance of deconvolution methods for various cell type removal scenarios compared to baseline with no cell type missing. Cell type reassignment metrics calculates assignment of missing celltype proportions from reference dataset.
>>>>>>> parent of c830ac4 (Updated scripts for better reproducibility)

**5\_Post\_analysis\_work**: The directory comprises scripts to generate the figures in the manuscript.

**6\_Quick\_Results**: The directory comprises scripts to generate results for the selective methods and analysis for quicker reproducibility.

<br>

#### How can you reproduce the results in the manuscript?

1. Software environment for R and Python used during the original analysis can be found under the`/0_SoftwareEnvironment` directory.


2. For reproducibility, the project utilises renv functionality. To setup renv functionalities for your analysis, please follow the steps below;
	
	i. Read [getting started with renv](https://rstudio.github.io/renv/articles/renv.html) carefully to understand the functionalities provided by the renv R package.
	
	ii. Creating an R project in RStudio: `RStudio > File > New Project > Existing Directory > Select the /Processing` directory in the sFSS. You can refer to [RStudio projects](https://support.posit.co/hc/en-us/articles/200526207-Using-RStudio-Projects) for more details on how the projects in RStudio work.
	
<<<<<<< HEAD
	**III**. Executing `Renv_Setup.R` under the `/Processing` directory to initialise the renv infrastructure for the R project:
	- The `Renv_Setup.R` script expects a **GITHUB_PAT token (classic)** to be set in the R environment. Please make sure you edit the script to add your own token. More details on the github token can be found [here](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens).
	- If using RStudio: execute the script line-by-line, read the comments within the scripts.
	- If using bash/shell terminal: <br>
		&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue" size=4>`Rscript Renv_setup.R`</font> 
=======
	iii. Open `Renv_Setup.R` under the `/Processing` directory in the newly created R project in RStudio and execute the script line-by-line.
>>>>>>> parent of c830ac4 (Updated scripts for better reproducibility)
		
		The required files (renv.lock, .Rprofile, renv/activate.R) should already be in your R project directory;
		If not, ensure it resides in the same directory as the Renv_setup.R file.
	
<<<<<<< HEAD
	> **<font color="red">Note:</font>** A collaborator can still experience minor discrepancies in the results while using the renv functionality due to the platform(OS) dependencies.
	
	<br>
	
	<font color="darkblue" size=4>§ Steps to set up conda functionalities for your analysis </font> <br>
	
	**I**. 	<font color="brown">[optional]</font> Read [<font color="blue">getting started with conda</font>](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html) carefully to understand the functionalities provided by the conda environment.
	
	**II.**	If you do not already have conda installed on your machine, please install it by following [<font color="blue"> steps to install conda</font>](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).
	
	**III**. To create a conda virtual environment similar to the analysis use the <font color="brown">*environment.yml*</font> file available in the root directory or under `0_SoftwareEnvironment/Python` directory. The command for creating conda virtual environment from command line terminal/prompt is as below;
	
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue" size=4>`conda env create -f environment.yml`</font> 

	> **<font color="red">Note:</font>** Depending on your platform(OS) installation of few python packages/libraries might 	fail, please check the logs/status of the environment creation carefully.
	
<br>
<br>

#### <font color="darkblue">Module-1: Navigate to "1\_Generate\_sc\_ref_data/Code/" directory </font>

* Generating a single-cell reference dataset from multiple single-cell datasets, UMAP representations included in the supplementary information of the manuscript using the commands as below; 
	
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Init_env.R`</font> <br>
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript SC_ref_data.R` </font> 
	
* Generating single-cell reference datasets for all removal scenarios based on cell type removal using the commands as below;
	
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Init_env.R`</font> <br>
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript SC_ref_data_scenarios.R`</font>
	> <font color="red">Note:</font> This script requires the single-cell reference dataset generated in the previous step.
	
<br>

#### <font color="darkblue">Module-2: Navigate to "2\_Simulating\_ST\_data/Code/" directory </font>
=======
	> **Note**: a collaborator can still experience minor discrepancies in the results while using the renv functionality due to the package dependencies, which are unavailable or non-functional and had to be updated to a newer version.
	
	> The differences in the package versions can be cross checked by the session info file provided under the `/Settings` directory in the every sub-analysis directory (e.g., `/1_Generate_sc_ref_data`) with the output of `sessionInfo()` command.


>>>>>>> parent of c830ac4 (Updated scripts for better reproducibility)

3. To load the required R packages in the session and set up the path to results, data and R/Python script directories.
	
	*R script:* `1_Generate_sc_ref_data/Code/Init_env.R`

<<<<<<< HEAD
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Init_env.R`</font> <br>
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Generate_ST_data.R`</font>
	
	> &emsp;&emsp;&emsp;&emsp;<font color="red">Note:</font> R script expects five command line arguments in the order mentioned below;<br>
	>	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;a. min number of cell types per spot <br>
	>	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;b. max number of cell types per spot <br>
	>	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;c. min number of cells per spot <br>
	>	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;d. max number of cells per spot <br>
	>	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;e. index of simulated ST dataset <font color="brown", size=2>[options: 1, 2, 3] </font> <br>
=======
4. To generate a single-cell reference dataset from multiple single-cell datasets. The paths to the processed data folder are already set through the 'Init_env.R' script.

	*R script:* `1_Generate_sc_ref_data/Code/SC_ref_data_donorwise.R`
	> *Note: This script takes a long time when executed on a local machine with (16GB RAM). This step can be skipped using the final result at `1_Generate_sc_ref_data/Results/sc.ref.data.rds`.*
If you do not have the `sc.ref.data.rds` file, please contact the developer to get it.

5. To generate single-cell reference datasets for different scenarios based on removing cell types.
	
	*R script:* `1_Generate_sc_ref_data/Code/SC_ref_data_scenarios.R`
	> *Note: This script requires the single-cell reference dataset generated in the previous step.*

6. To generate simulated spatial transcriptomics datasets using single-cell reference data. The four datasets vary in terms of the number of cells and cell types present per spatial location.

	*R script:* `2_Simulating_ST_data/Code/Generate_ST_data.R`
	> **Note: R script expects 5 command line arguments in the order mentioned below;**  
	> 	*a. min number of cell types per spot <br>
	> 	b. max number of cell types per spot <br>
	>	c. min number of cells per spot <br>
	>	d. max number of cells per spot <br>
	>	e. index of simulated ST dataset (either 1/2/3) <br>*
>>>>>>> parent of c830ac4 (Updated scripts for better reproducibility)
	
		E.g., Commands to generate 1st, 2nd & 3rd simulated ST datasets,
		Rscript Generate_ST_data.R 4 8 10 15 1
		Rscript Generate_ST_data.R 1 5 10 15 2
		Rscript Generate_ST_data.R 1 5 3 7 3

<<<<<<< HEAD
<br>

#### <font color="darkblue">Module-3: Navigate to "3\_ST\_methods/Code/" directory </font>
=======
7. To predict cell type proportions using the cell type deconvolution methods. The shell scripts execute all deconvolution methods simultaneously in parallel or serial combination. The standalone scripts for each deconvolution method also exist and expects the same arguments as that of shell script.
>>>>>>> parent of c830ac4 (Updated scripts for better reproducibility)

	*Shell scripts to execute R-based deconvolution methods:* `3_ST_methods/Code/Execute_R_based_methods.sh`

<<<<<<< HEAD
	Shell/batch scripts to execute R/Python-based deconvolution methods: <br>
	<font color="blue">`Execute_R_based_methods.sh`</font>, <br>
	<font color="blue">`Execute_python_based_methods.sh`</font>, <br>
		<font color="blue">`Execute_R_based_methods.bat`</font>, <br>
	<font color="blue">`Execute_python_based_methods.bat`</font>

	> <font color="red">Note:</font> The script to execute R-based/Python-based methods executes each method in parallel for a provided number of reference datasets and removal scenarios; the required command line arguments are as below:<br>
	> &emsp;&emsp;&emsp;&emsp;a. index of the ST dataset <font color="brown", size=2>[options: 1, 2, 3] </font> <br>
	> &emsp;&emsp;&emsp;&emsp;b. total number of single-cell reference datasets <font color="brown", size=2>[options per removal scenario; <br>
	&emsp;&emsp;&emsp;&emsp;&emsp;rm0= 1:1, rm1= 1:13, rm2= 1:5, rm3= 1:5, rm5= 1:1, rm10= 1:5, rm11= 1:5] </font> <br>
	> &emsp;&emsp;&emsp;&emsp;c. removal scenario <font color="brown", size=2>[options: rm0, rm1, rm2, rm3, rm5, rm10, rm11] </font>
		
		e.g. command line input for removal of no cell type removal (baseline) & one or more cell type removal scenarios would be then as follows,
		(1) ./Execute_R_based_methods.sh 1 1 rm0
		(2) ./Execute_R_based_methods.sh 1 13 rm1
		(3) ./Execute_R_based_methods.sh 1 5 rm2
		(4) ./Execute_R_based_methods.sh 1 5 rm3
		(5) ./Execute_R_based_methods.sh 1 1 rm5
		(6) ./Execute_R_based_methods.sh 1 5 rm10
		(7) ./Execute_R_based_methods.sh 1 5 rm11
		
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="brown", size=2> replace `Execute_R_based_methods.sh` by `Execute_python_based_methods.sh` for executing python-based methods. </font>
	
	*The standalone scripts for each deconvolution method expect the command line arguments as below; notice the different versions of argument <font color="red">b </font>*.

	> *a. index of the ST dataset <br>
	> b. <font color="blue", size=2>for R-based methods: </font> index of single-cell reference datasets [varies per removal scenario] <br>
	&emsp;&emsp;&emsp;&emsp;<font color="blue", size=2>for Python-based methods: </font> total number of single-cell reference datasets <br>
	> c. path to simulated ST datasets directory <br>
	> d. path to single-cell reference datasets directory <br>
	> e. path to save deconvolution results directory* <br>
=======
	> **Note: The shell script executes R-based methods in parallel and in a loop <br>**
	> *a. total number of ST datasets [options: 1, 2, 3] <br>
	> b. total number of single-cell reference datasets [options varies per removal scenario; <br> rm0= 1, rm1= 13, rm2= 5, rm3= 5, rm5= 1, rm10= 5, rm11= 5] <br>
	> c. path to simulated ST datasets directory <br>
	> d. path to single-cell reference datasets directory <br>
	> e. path to save deconvolution results directory <br>*
		
		E.g. commands for removal of no cell type mismatch (baseline) & one or more cell type mismatch scenarios would be then as follows,
		
		(1) ./Execute_R_based_methods.sh 1 1 "../../2_Simulating_ST_data/Results/Spatial_Data/" "../../1_Generate_sc_ref_data/Results/" "../Results/rm0/"
		(2) ./Execute_R_based_methods.sh 1 13 "../../2_Simulating_ST_data/Results/Spatial_Data/" "../../1_Generate_sc_ref_data/Results/rm1/" "../Results/rm1/"
		(3) ./Execute_R_based_methods.sh 1 5 "../../2_Simulating_ST_data/Results/Spatial_Data/" "../../1_Generate_sc_ref_data/Results/rm2/" "../Results/rm2/"
		(4) ./Execute_R_based_methods.sh 1 5 "../../2_Simulating_ST_data/Results/Spatial_Data/" "../../1_Generate_sc_ref_data/Results/rm3/" "../Results/rm3/"
		(5) ./Execute_R_based_methods.sh 1 1 "../../2_Simulating_ST_data/Results/Spatial_Data/" "../../1_Generate_sc_ref_data/Results/rm5/" "../Results/rm5/"
		(6) ./Execute_R_based_methods.sh 1 5 "../../2_Simulating_ST_data/Results/Spatial_Data/" "../../1_Generate_sc_ref_data/Results/rm10/" "../Results/rm10/"
		(7) ./Execute_R_based_methods.sh 1 5 "../../2_Simulating_ST_data/Results/Spatial_Data/" "../../1_Generate_sc_ref_data/Results/rm11/" "../Results/rm11/"
		

	*Shell scripts to execute python-based deconvolution methods:* `3_ST_methods/Code/Execute_python_based_methods.sh`
>>>>>>> parent of c830ac4 (Updated scripts for better reproducibility)
	
	> **Note: The shell script executes Python-based methods in parallel <br>**
	> *a. index of the ST dataset [options: 1, 2, 3] <br>
	> b. total number of single-cell reference datasets [options varies per removal scenario; <br> rm0= 1, rm1= 13, rm2= 5, rm3= 5, rm5= 1, rm10= 5, rm11= 5] <br>
	> c. path to simulated ST datasets directory <br>
	> d. path to single-cell reference datasets directory <br>
	> e. path to save deconvolution results directory <br>*

	*The standalone scripts for each deconvolution method also expect same command line arguments as above*<br>

<<<<<<< HEAD
<br>

#### <font color="darkblue">Module-4: Navigate to "4\_Analysis\_results/Code/" directory </font>

* Estimating performance metrics for all the removal of cell types scenarios.

	> <font color="red">Note:</font> Estimating performance metrics needs all the deconvolution results generated in Module-3.

    - **For JSD calculations**:
   
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Init_env.R`</font> <br> 
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Get_JSD_results.R `</font>
	> &emsp;&emsp;&emsp;&emsp;<font color="red">Note: </font> The required command line arguments are as below:<br>
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;a. index of the simulated ST dataset <font color="brown", size=2>[options: 1, 2, 3] </font> <br>	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;b. <font color="brown", size=2>[optional; default is all scenarios] </font> removal scenarios separated by comma (<font color="brown", size=2>rm0 </font> is included by default) <font color="brown", size=2>[options: rm1, rm2, rm3, rm5, rm10, rm11] </font> <br>
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;c. <font color="brown", size=2>[optional; default is all methods] </font> name of the deconvolution methods separated by comma <br>
		
		e.g.1, command line input for calculating JSD metrics for methods CARD and RCTD with 1st simulated ST data and removal of one, two & three cell type scenarios would be then as follows,
		Rscript Get_JSD_results.R 1 "rm1","rm2","rm3"  "CARD","RCTD"
		
		e.g.2, command line input for calculating JSD metrics for all methods with 1st simulated ST data and all removal of cell type scenarios would be then as follows,
		Rscript Get_JSD_results.R 1

   - **For RMSE calculations**:

   &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Init_env.R` </font> <br> 
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Get_RMSE_results.R` </font>
   > &emsp;&emsp;&emsp;&emsp;<font color="red">Note:</font> The required command line arguments are as below:<br>
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;a. index of the simulated ST dataset <font color="brown", size=2>[options: 1, 2, 3] </font> <br>	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;b. removal scenarios separated by comma (<font color="brown", size=2>rm0 </font> is included by default) <font color="brown", size=2>[options: rm1, rm2, rm3] </font> <br>
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;c. <font color="brown", size=2>[optional; default is all methods] </font> name of the deconvolution methods separated by comma <br>
		
		e.g., command line input for calculating RMSE metrics for methods CARD and RCTD with 1st simulated ST data and removal of one, two & three cell type scenarios would be then as follows,
		Rscript Get_RMSE_results.R 1 "rm1","rm2","rm3"

   - **For cell type reassignment calculations:** <br>
   
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Init_env.R`</font> <br> 
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Get_celltype_assignment_results.R `</font>
	> &emsp;&emsp;&emsp;&emsp;<font color="red">Note: </font> The required command line arguments are as below:<br>
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;a. index of the simulated ST dataset <font color="brown", size=2>[options: 1, 2, 3] </font> <br>	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;b. <font color="brown", size=2>[optional; default is all scenarios] </font> removal scenarios separated by comma (<font color="brown", size=2>rm0 </font> is included by default) <font color="brown", size=2>[options: rm1, rm2, rm3, rm5, rm10, rm11] </font> <br>
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;c. <font color="brown", size=2>[optional; default is all methods] </font> name of the deconvolution methods separated by comma <br>
		
		e.g., command line input for calculating cell type reassignment for method CARD and RCTD with 1st simulated ST data and removal of one, two & three cell type scenarios would be then as follows,
		Rscript Get_celltype_assignment_results.R "CARD","RCTD" 1 "rm1","rm2","rm3"    

<br>

#### <font color="darkblue">Module-5: Navigate to "5\_Post\_analysis\_work/Code/" directory </font>

* Generating figures/plots included in the manuscript for calculated JSD, RMSE and cell type reassignment estimates for the specified cell type removal scenarios.

	> <font color="red">Note:</font> this module expects calculations results from Module-4
	
   - Generate plots for cell type correlation and cell type reassignment plots for cell type removal scenarios rm1, rm2, and rm3 using commands as below;
   
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Init_env.R`</font> <br> 
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Celltype_correlation.R `</font> <br>
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Fig_celltype_assignments.R`</font> <br>
   	> &emsp;&emsp;&emsp;&emsp;<font color="red">Note:</font> <font color="blue">`Fig_celltype_assignments.R`</font> expects command line arguments as below:<br>
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;a. index of the simulated ST dataset <font color="brown", size=2>[options: 1, 2, 3] </font> <br>
   
   - Generate plots for JSD/RMSE plots for all cell type removal scenarios using commands as below;
   
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Init_env.R`</font> <br> 
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Fig_jsd_rmse_plots.R`</font> <br>
	> &emsp;&emsp;&emsp;&emsp;<font color="red">Note:</font> <font color="blue">`Fig_jsd_rmse_plots.R`</font> expects command line arguments as below; (<font color="brown" size=2>list of methods and removal scenarios are adapted from JSD/RMSE calculations in Module-4 </font>)<br>
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;a. index of the simulated ST dataset <font color="brown", size=2>[options: 1, 2, 3] </font> <br>

   - Generate summary funky plot for an overview of ranking of deconvolution methods across all the removal scenarios using commands as below;
   
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Init_env.R`</font> <br> 
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<font color="blue">`Rscript Funky_plots.R`</font> <br>
	> &emsp;&emsp;&emsp;&emsp;<font color="red">Note:</font> <font color="blue">`Funky_plots.R`</font> expects command line arguments as below;<br>
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;a. index of the simulated ST dataset <font color="brown", size=2>[options: 1, 2, 3] </font> <br>
=======
	

8. Estimate performance metrics for all the scenarios (baseline and removal) in `4_Analysis_results/Code`
    - *For JSD calculations:* `Get_JSD_results.R`
    - *For RMSE calculations:* `Get_RMSE_results.R`
    - *For cell tyepe reassignment calculations:* `Get_celltype_assignment_results.R`
    
	> *Note: Estimating performance metrics needs deconvolution results generated in the previous step.*

9. Generate figures and results for calculated JSD, RMSE and cell type reassignment estimates for all the cell type removal scenarios.

	*R scripts* under `5_Post_analysis_work/Code`
	
	`Fig_jsd_rmse_plots.R`, `Funky_plots.R`, `Plot_CT_assignment_rm1.R`, `Plot_CT_assignment_rm1.R`, `Plot_CT_assignment_rm1.R`, `Plot_CT_assignment_rm1.R`

>>>>>>> parent of c830ac4 (Updated scripts for better reproducibility)

<br>
<br>
### How to execute the 6\_Quick_Results module?
TBA


### Issues:

**1. Installation of 'fields' R package on macOS**

The analysis is carried out with fields package version 13.3 developed under R version 4.1 and RStudio built for x86\_64 architecture.

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

to resolve this error, 

- you can install the latest version 15.2 ([R-binaries](https://cran.r-project.org/web/packages/fields/index.html) available for arm\_64); this will affect the CARD package since it needs to be updated to the latest version as well.

- install R and RStudio built for x86\_64 architecture and run the analysis using it. If you still get the same error, the compiler uses the default arm\_64 architecture to install R packages that need compilation.

<br>

**2. Compilation failed error on Windows OS**

R and RStudio on Windows OS require rtools (a toolchain for building R and R packages). The correct version of rtools can be downloaded from [here](https://cran.r-project.org/bin/windows/Rtools/).

	Ensure you do not have white spaces in your paths while installing the rtools.
	
<br>

**3. Miniconda installation prompt**

		No non-system installation of Python could be found.
		Would you like to download and install Miniconda?
		Miniconda is an open-source environment management system for Python.
		See https://docs.conda.io/en/latest/miniconda.html for more details.
		 
		Would you like to install Miniconda? [Y/n]:

If you see above prompt in your execution of analysis, press `n` and stop the execution.

Before moving forward please check your python installation before moving forward with `Sys.which("python")` or `Sys.which("python3")` in RStudio console.

- If you see `" "`, please install the python version 3.9.7 and try again.
- If you see `"/path_to_python_installation"`, please contact the developer or create an issue on the GitHub repository with the error message, R session info, and python installation details (path, version, etc).

		
