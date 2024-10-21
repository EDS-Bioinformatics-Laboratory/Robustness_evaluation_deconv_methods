# Analysing the deconvolution results



*<u>**Instructions**</u>*

* *Remove instructions once you finalized the content of this file.*
* *Subdirectories*
  * ***Code**. Contains the (in-house developed) software for the computational analysis.*
  * ***CodeDocumentation**. Any documentation that is available about the software such as requirement specifications, software design, source code documentation, testing requirements, and end-user instructions. Alternatively, part of this documentation can be placed in the code itself, or in the code 0_README.md file.*
  * ***Data**. Raw, meta, and pre-processed data.*
  * ***Notebooks**. Contains interactive notebooks (e.g, R or Jupyter notebooks), which may include data and results. Minimize the size of the Notebook prior to pushing to GitHub.*
  * ***Results**. (Intermediate) results such as figures and tables from the computational analysis. Additional sub-directories within /Results, to organize the results, are allowed.* 
  * ***Settings**. Contains any file with input parameters for e.g., statistical analysis and computational simulations. If only few parameters or parameter files are required then these might also be placed in the \Code directory or directly in the code itself.*



*<u>Manual steps:</u>*

***Avoid** manual data manipulation steps. Instead of manually changing data (e.g., format conversion, filtering), aim to make use of small custom scripts. Manual steps are a main cause for irreproducible results.*  

* *Having an executable description of the sequence of steps taken helps to reproduce results (e.g., Snakemake (https://snakemake.github.io/) or pytask (https://github.com/pytask-dev/pytask)).*



*<u>Intermediate results and hierarchical output</u>*

*Record intermediate results (preferable in a standardized format). Generate hierarchical analysis output, allowing layers of increasing detail to be inspected. This can reveal discrepancies toward what is assumed, and can in this way uncover bugs or faulty interpretation that are not apparent in the final results. It also allows any inconsistency to be tracked to the step where the problem occurs.  It also allows critical examination of the full process behind a result.*



*<u>Regeneration of figures and tables:</u>* 

*For any figure or table that ends up in a publication, report, or presentation at meeting, the underlying data and a stand-alone piece of code should be available to regenerate the figure. It also allows easy modification of a figure and to retrieve the data of the figure (instead of having to redo a complete analysis). Equally important, the data of the figure can be further analysed or inspected.*



*== END INSTRUCTIONS ==*



**Operating System(s) / version(s) used during development (and testing):**

MacOS Ventura 13.4.1


**Software environment:** 

- R version : 4.1.2 (2021-11-01)
- Rstudio IDE: RStudio 2023.03.0+386 "Cherry Blossom"



####**Conceptual description of methodology**

We have used Jansen-Shannon divergence (JSD) and the root-mean-square error (RMSE) as the two similarity metrics to understand the performance of robustness of the deconvolution methods. We evaluated if a deconvolution method's predicted cell type proportion accurately represents true proportion as per the ground truth in the three simulated spatial transcriptomics data. Along with these two-similarity metrics, we have also calculated the cell type assignment, to understand how the proportion of the missing cell type is handled by the deconvolution methods.

**Jensen-Shannon divergence (JSD)**: JSD describes the similarity between two probability distributions (here, ground truth and the predicted proportion in a spot). The JSD is a symmetrized and smoother version of the Kullback-Liebler divergence (KLD). We used the following equation to calculate the JSD value for each spot.

$$
KLD\ (pred\ |\ act) = \sum_{x\ ∈\ X} pred(x)\ *\ log\left({pred(x)\over{act(x)}} \right)
$$
$$
M\ =\ {{1}\over{2}}\ {(pred\ + act)}
$$
$$
JSD\ (pred\ |\ act)\ =\ {{1}\over{2}}\ {KLD\ (pred\ |\ M)}\ +\ {{1}\over{2}}\ {KLD\ (act\ |\ M)}
$$
	here, pred and act are the predicted and ground truth proportion of cell types in a spot respectively defined in probability space X. If the two distributions are identical, the JSD value is zero. The lower the value of JSD, the better the predictions of cell type deconvolution are.

**Root Mean Square Error (RMSE)**: RMSE represents the measurement of the difference between a predicted proportion and the ground truth proportion of cell type in a spot and not between the spots. We used the following equation to calculate the RMSE for each spot.

$$
RMSE = \sqrt{\sum_{i = 1}^{n}(predicted_i - actual_i)\over{n}}
$$
	here, predicted<sub>i</sub> and actual<sub>i</sub> are the predicted and ground truth proportion for cell type i in the spot respectively and n is the number of distinct cell types present in the data. The lower the value of RMSE, the better the predictions of cell type deconvolution are.

**Cell type reassignment**: Cell type reassignment is calculated using the change in baseline proportions of all cell type present in the reference data. We considered the proportions of cell types across all the spots in ST data when no cell types are missing in the reference data (i.e., baseline) and when one or more cell types are missing. The difference between these proportions has been normalized to one to estimate the proportion reassignment. The procedure is as mentioned below.

		Let B be the predicted cell type proportions at the baseline scenario;
		Let P be the predicted cell type proportions at the removal scenario;
		
		For all spots in the ST data
		    Let X be the change in proportions;
		    Let c be the common cell types in both the predicted proportions;
		    Let m be the missing cell type(s) in the predicted proportion for removal scenario;
	
		    For all c
		        X = P – B
		    Normalize X for proportion of m for baseline scenario


