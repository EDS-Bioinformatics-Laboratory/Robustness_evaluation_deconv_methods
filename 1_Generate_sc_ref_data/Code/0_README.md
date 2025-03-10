### CODE


**Programming language + version:**

- R version 4.1.2 (2021-11-01)
- Other packages used during the R session can be found in `sessionInfo.txt` file under `/Settings` directory


**What types of documentation did you produce:**

Each script produces some intermediate/end results and can be found in Results directory in this analysis.

**Which integrated development environment (IDE) did you use:**

RStudio 2023.03.0+386

**How did you test the code:**

For manual testing, each script can be run independently. More instructions are in another README.md and are available in the project's `/Processing` directory. The scripts should be executed in the order mentioned below.

**Init_env.R:** (optional if running in terminal or command prompt) The script loads the required R libraries (packages) for the analysis. It also write the session info collecting all the R package versions used during analysis along with collection of citations of the packages. 


**SC\_ref\_data.R:** The Rscript reads two scRNA-seq datasets (lymph node stromal cells dataset and tabula sapiens consortium's lymph node dataset). The tabula sapiens lymph node dataset is split donor-wise into three individual datasets out of which one dataset is further split platform-wise. The script further takes care quality control check, pre-processing, integration of datasets, annotations of the integrated dataset, and downsampling of the integrated data to generate a single cell reference dataset.

**SC\_ref\_data_scenarios.R:** It generates single cell reference data versions by removing one or more cell types from the earlier reference dataset via `SC_ref_data_donorwise.R` script. (should be executed after the scRNA-seq reference dataset is generated) 

**logs:** This consists of log (text) files from the execution for reference.