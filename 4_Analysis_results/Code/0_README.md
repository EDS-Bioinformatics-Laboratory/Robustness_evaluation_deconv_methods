### CODE


**Programming language + version:**

R version 4.1.2 (2021-11-01)

**What types of documentation did you produce:**

Each script produces some intermediate/end results and can be found in Results directory in this analysis.

**Which integrated development environment (IDE) did you use:**

RStudio 2023.03.0+386

**How did you test the code:**

Manual testing, each script can be run independently. More instructions can be found in another README.md and can be found in root folder of this analysis.

**Init_env.R:** The script loads the required R libraries (packages) for the analysis. It also write the session info collecting all the R package versions used during analysis along with collection of citations of the packages. 

**Get\_celltype\_assignment_results.R:** The missing cell type proportion is reassigned to available cell types in the reference data. This script estimates the cell type reassignment.

**Get\_JSD\_results.R:** Calculates the difference in JSD values for all cell type removal scenarios compare to baseline scenarion with no missing cell type. The script calculates mean JSD values across a spatial location.

**Get\_RMSE\_results.R:** Calculates the difference in RMSE values for all cell type removal scenarios compare to baseline scenarion with no missing cell type. The script calculates mean RMSE values across a spatial location.


> <font color="red">Note:</font> More details for every script are available in the header of the script.