### CODE


**Programming language + version:**

R version 4.1.2 (2021-11-01)

**What types of documentation did you produce:**

Each script produces some intermediate/end results and can be found in Results directory in this analysis.

**Which integrated development environment (IDE) did you use:**

RStudio 2023.03.0+386

**How did you test the code:**

Manual testing, each script can be run independently. More instructions can be found in another README.md and can be found in root folder of this analysis.

**Init_env.R:** The script loads the required R libraries (packages) for the analysis. It also writes the session info, collecting all the R package versions used during analysis.

**Generate\_ST_data.R:** Generate simulated ST datasets using the scRNA-seq dataset as a reference dataset. The script requires command line arguments as mentioned below:

1. min number of cell types per spot
2. max number of cell types per spot
3. min number of cells per spot
4. max number of cells per spot
5. simulated ST dataset number (either 1/2/3)
