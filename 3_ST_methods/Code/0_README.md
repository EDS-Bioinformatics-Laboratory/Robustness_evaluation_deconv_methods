### CODE



<u>***Instructions***</u>

* *Remove instructions once you finalized the content of this file.*



***Clean code*** 

*Write programs for people, not computers (code should be easily read and understood). Adopting clean code practices helps to standardize and organize software code in order to enhance readability. This allows developers to concentrate on core functionality and reduce errors.* 

*A few general rules:*

- *Reduce complexity. Decompose programs into functions and modules (instead of making very long scripts).*
- *Limit number of function arguments and length of functions.*
- *Be ruthless about eliminating duplication (use functions).*
- *Always search for well-maintained software libraries that do what you need, and test these libraries     before relying on them.*
- *Do not comment and  uncomment sections of code to control a program's behaviour. This is a recipe for making your software irreproducible.* 
- *Choose a style guide—And stick with it (see for example, https://peps.python.org/pep-0008/ or https://web.stanford.edu/class/cs109l/unrestricted/resources/google-style.html)*
- *Give functions and variables meaningful names.*
- *Make dependencies and requirements explicit*
- *Refactor as needed (process of restructuring your code without changing its interface).*

- *Write error messages that provide solutions or point to your documentation.*



***Code testing***

*Interrogate specific and isolated coding behaviour to reduce coding errors and ensure intended functionality, especially as code increases in complexity. Describe if software tests performed and how to re-run these tests.*



***Integrated development environment** (IDE)*

*Consider using an IDE (integrated development environment (IDE; )e.g., PyCharm and DataSpell from JetBrains (https://www.jetbrains.com/), RStudio (https://www.rstudio.com/)* 



***Software documentation***

*Not all types of software documentation on this list need to be written for a single research or support project but a combination of several should be selected depending on the nature of the project. Whether more or less documentation is needed for your project will depend on its scale and complexity. For research project we need at least the (external) source code documentation and user documentation.*

- *Requirements Specification*

- *Software Design*

- ***(External) source code documentation***

  - *Place a brief explanatory comment at the start of every program.*

  - *Document the input and output of your script/functions.*
  - *Describe what the code is doing and why.*

- *Testing Requirements*
- ***End-User Instructions (including a quick-start guide)***
  - *How to install and configure the software.*
  - *Where to find its full documentation.*
  - *Include examples how to execute the code to produce certain output. Provide a simple example or test dataset.*
  - *Under what license it’s released.*
- *Include a help command for command line interfaces*



***Write comments as you code (not afterwards).** Modern integrated development environments (IDEs) will often automatically generate documentation strings as you write code, which removes the burden of having to remember to write comments. (e.g., PyCharm, DataSpell).*



***Version control your documentation.** Keep your documentation inside your Git repository along with the rest of your files. In addition/alternatively ‘Read the Docs’ (https://readthedocs.org/) provides an approach for hosting and automatically building documentation.* 



***Use automated documentation tools***

*For example, using Python docstrings (https://realpython.com/documenting-python-code/) together with Sphinx https://www.sphinx-doc.or) to automatically generate documentation.* 



 *== END INSTRUCTIONS ==*



**Programming language + version:**

R version 4.1.2 (2021-11-01)

**What types of documentation did you produce:**

Each script produces some intermediate/end results and can be found in Results directory in this analysis.

**Which integrated development environment (IDE) did you use:**

RStudio 2023.03.0+386

**How did you test the code:**

Manual testing, each script can be run independently. More instructions can be found in another README.md and can be found in root folder of this analysis.

**Init_env.R:** The script loads the required R libraries (packages) for the analysis. It also write the session info collecting all the R package versions used during analysis along with collection of citations of the packages. 


**Execute\_python\_based\_methods.sh:** Script executes python based deconvolution methods (Cell2Location and Stereoscope) in serial combination at once for the provided user inputs as the command line arguments.


**Execute\_R\_based\_methods.sh:** Script executes python based deconvolution methods (CARD, MuSiC, RCTD, SCDC, Seurat and SPOTlight) in parallel combination at once for the provided user inputs as the command line arguments.

**CARD.R:** Standalone R script for execution of CARD, a deconvolution method tailored for spatial transcriptomics data. Requires command line arguments.

**Cell2Location.py:** Standalone python script for execution of Cell2Location, a deconvolution method tailored for spatial transcriptomics data. Requires command line arguments.

**MuSiC.R:** Standalone R script for execution of MuSiC, a deconvolution method tailored for bulk RNA-sequencing data. Requires command line arguments.

**RCTD.R:** Standalone R script for execution of RCTD, a deconvolution method tailored for spatial transcriptomics data. Requires command line arguments.

**SCDC.R:** Standalone R script for execution of SCDC, a deconvolution method tailored for bulk RNA-sequencing data. Requires command line arguments.

**Seurat.R:** Standalone R script for execution of Seurat, a deconvolution method tailored for spatial transcriptomics data. Requires command line arguments.

**SPOTlight.R:** Standalone R script for execution of SPOTlight, a deconvolution method tailored for spatial transcriptomics data. Requires command line arguments.

**Stereoscope.R:** Standalone python script for execution of Stereoscope, a deconvolution method tailored for spatial transcriptomics data. Requires command line arguments.


> *Note: All scripts (except Init_env.R) expects 5 command line arguments in order mentioned below;*<br>

	a. total number of ST datasets [options: 1, 2, 3, 4],
	b. total number of single-cell reference datasets [options varies per removal scenario; rm1= 13, rm2= 5, rm3= 5, rm5= 1, rm10= 5, rm11= 5],
	c. path to simulated ST datasets directory,
	d. path to single-cell reference datasets directory,
	e. path to save deconvolution results directory
	
	E.g. script for removal of one cell type scenario,
	./Execute_R_based_methods.sh 1 13 "../../2_Simulating_ST_data/Results/Spatial_Data/" "../../1_Generate_sc_ref_data/Results/rm1/" "../Results/rm1/"
	Rscript CARD.R 1 13 "../../2_Simulating_ST_data/Results/Spatial_Data/" "../../1_Generate_sc_ref_data/Results/rm1/" "../Results/rm1/"