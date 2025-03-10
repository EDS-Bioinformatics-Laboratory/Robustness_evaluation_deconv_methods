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


**Execute\_python\_based\_methods.sh/.bat:** Script executes python based deconvolution methods (Cell2Location and Stereoscope) in serial combination at once for the provided user inputs as the command line arguments.


**Execute\_R\_based\_methods.sh/.bat:** Script executes python based deconvolution methods (CARD, MuSiC, RCTD, SCDC, Seurat and SPOTlight) in parallel combination at once for the provided user inputs as the command line arguments.

**CARD.R:** Standalone R script for execution of CARD, a deconvolution method tailored for spatial transcriptomics data. Requires command line arguments.

**Cell2Location.py:** Standalone python script for execution of Cell2Location, a deconvolution method tailored for spatial transcriptomics data. Requires command line arguments.

**MuSiC.R:** Standalone R script for execution of MuSiC, a deconvolution method tailored for bulk RNA-sequencing data. Requires command line arguments.

**RCTD.R:** Standalone R script for execution of RCTD, a deconvolution method tailored for spatial transcriptomics data. Requires command line arguments.

**SCDC.R:** Standalone R script for execution of SCDC, a deconvolution method tailored for bulk RNA-sequencing data. Requires command line arguments.

**Seurat.R:** Standalone R script for execution of Seurat, a deconvolution method tailored for spatial transcriptomics data. Requires command line arguments.

**SPOTlight.R:** Standalone R script for execution of SPOTlight, a deconvolution method tailored for spatial transcriptomics data. Requires command line arguments.

**Stereoscope.R:** Standalone python script for execution of Stereoscope, a deconvolution method tailored for spatial transcriptomics data. Requires command line arguments.

**logs:** This consists of log (text) files from the execution for reference.


> *<font color="red">Note</font>: A detail explanation about how to execute these scripts is provided in the 0_README.md file available in the root (/Processing) directory*<br>
