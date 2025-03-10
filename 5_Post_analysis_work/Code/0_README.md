### CODE


**Programming language + version:**

R version 4.1.2 (2021-11-01)

**What types of documentation did you produce:**

Each script produces some intermediate/end results and can be found in Results directory in this analysis.

**Which integrated development environment (IDE) did you use:**

RStudio 2023.03.0+386

**How did you test the code:**

Manual testing, each script can be run independently. More instructions can be found in another README.md and can be found in root folder of this analysis.

**Init_env.R:** The script loads the required R libraries (packages) for the analysis. It also write the session info collecting all the R package versions used during analysis.

**Celltype_correlation.R:** The script calculates and plots the cell type ranking based on correlation between the gene expression profilesof the cell types present in the single cell reference data.


**Fig\_celltype_assignments.R:** The script generats cell type reassignment plots for cell type removal scenarios rm1, rm2, and rm3. More details are within the script in the header.

**Fig\_jsd\_rmse_plots.R:** The script generates plots for JSD/RMSE plots for all cell type removal scenarios. More details are within the script in the header.

**Funky_plots.R:** The script generates summary funky plot for an overview of ranking of deconvolution methods across all the removal scenarios.

**Validate_numbers.R:** The script generates barplot for number of cells per celltype in the integrated single-cell reference dataset. It also generates heatmap plots for JSD/RMSE values for all scenarios in each ST dataset. These results can be cross-verified with the results produced in the original analysis.


**logs:** This consists of log (text) files from the execution for reference.