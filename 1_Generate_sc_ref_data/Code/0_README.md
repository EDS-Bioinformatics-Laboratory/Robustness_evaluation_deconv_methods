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

- R version 4.1.2 (2021-11-01)
- Other packages used during the R session can be found in `sessionInfo.md` file under `/Settings` directory


**What types of documentation did you produce:**

Each script produces some intermediate/end results and can be found in Results directory in this analysis.

**Which integrated development environment (IDE) did you use:**

RStudio 2023.03.0+386

**How did you test the code:**

For manual testing, each script can be run independently. More instructions are in another README.md and are available in the project's `/Processing` directory. The scripts should be executed in the order mentioned below.

**Init_env.R:** (optional if running in terminal or command prompt) The script loads the required R libraries (packages) for the analysis. It also write the session info collecting all the R package versions used during analysis along with collection of citations of the packages. 


**SC\_ref\_data_donorwise.R:** The Rscript reads two scRNA-seq datasets (lymph node stromal cells dataset and tabula sapiens consortium's lymph node dataset). The tabula sapiens lymph node dataset is split donor-wise into three individual datasets. The script further takes care quality control check, pre-processing, integration of datasets, annotations of the integrated dataset, and downsampling of the integrated data to generate a single cell reference dataset.

**SC\_ref\_data_scenarios.R:** It generates single cell reference data versions by removing one or more cell types from the earlier reference dataset via `SC_ref_data_donorwise.R` script. (should be executed after the scRNA-seq reference dataset is generated) 

**Clustree_map.R:** Verifies the cell type annotations of final scRNA-seq reference dataset generated via `SC_ref_data_donorwise.R` script using clustree maps and alluvial/sankey plots.

**SC\_ref\_data\_platform_wise.R:** This script generates necessary UMAP plots for supplementary information. Here, the tabula sapiens lymph node dataset is split platform-wise and integrated with the in-house LNSC dataset to understand the batch effects we observe when the dataset is split donor-wise.
