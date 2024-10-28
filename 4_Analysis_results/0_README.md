# Analysing the deconvolution results



**Operating System(s) / version(s) used during development (and testing):**

MacOS Ventura 13.4.1


**Software environment:** 

- R version : 4.1.2 (2021-11-01)
- Rstudio IDE: RStudio 2023.03.0+386 "Cherry Blossom"



####**Conceptual description of methodology**

We have used Jansen-Shannon divergence (JSD) and the root-mean-square error (RMSE) as the two similarity metrics to understand the performance of robustness of the deconvolution methods. We evaluated if a deconvolution method's predicted cell type proportion accurately represents true proportion as per the ground truth in the three simulated spatial transcriptomics data. Along with these two-similarity metrics, we have also calculated the cell type assignment, to understand how the proportion of the missing cell type is handled by the deconvolution methods.

**Jensen-Shannon divergence (JSD)**: JSD describes the similarity between two probability distributions (here, ground truth and the predicted proportion in a spot). The JSD is a symmetrized and smoother version of the Kullback-Liebler divergence (KLD). We used the following equation to calculate the JSD value for each spot.

$$
KLD\ (pred\ ||\ gt) = \sum_{k=1}^{C} pred(k)\ *\ log_2\left({pred(k)\over{gt(k)}} \right)
$$
$$
M\ =\ {{1}\over{2}}\ {(pred\ + gt)}
$$
$$
JSD\ (pred\ ||\ gt)\ =\ {{1}\over{2}}\ {KLD\ (pred\ ||\ M)}\ +\ {{1}\over{2}}\ {KLD\ (gt\ ||\ M)}
$$
	where *pred* and *gt* denote the predicted and ground truth cell type proportions of *C* cell types, respectively. The JSD is bounded between 0 and 1, and equal to zero if and only if the two distributions are identical. Therefore, a lower JSD indicates better performance.

**Root Mean Square Error (RMSE)**: RMSE represents the measurement of the difference between a predicted proportion and the ground truth proportion of cell type in a spot and not between the spots. We used the following equation to calculate the RMSE for each spot.

$$
RMSE = \sqrt{\sum_{k = 1}^{C}(pred_k - gt_k)^2\over{n}}
$$
	The lower RMSE indicates better performance.

**Cell type reassignment**: The predicted proportions of a cell type at baseline become zero if the cell type is removed from the reference data in one of the cell type mismatch scenarios. The decrease in the predicted proportions of the removed cell type(s) has to be accompanied by an increase in the predicted proportions of one or more of the remaining cell types. Therefore, we introduce the reassignment value to capture to which cell type(s) the proportions of the removed cell type(s) are reassigned. The calculation of the reassignment value is outlined below,


**B**: 1600 x 13 matrix of predicted cell type proportions at baseline for 1600 spots and 13 cell types

**P**: 1600 x (13 – *j*) matrix of predicted cell type proportions with *j* ∈ {1,2,3,5,10,11} cell types removed from the reference data

**S**: set of remaining cell types, that is |*S*| = 13 – *j*

**R**: set of removed cell types, that is |*R*| = *j*

**for** *s* **in** *S*

&emsp;&emsp;  _Calculate change in proportions w.r.t baseline for cell type **s** for all spots_
    
&emsp;&emsp;  **d** = **P**[:, *s*] - **B**[:, *s*]   

&emsp;&emsp;  _Calculate the reassignment value **rv** for cell type **s** when removing_

&emsp;&emsp;  _Here, the numerator corresponds to the total change in proportions w.r.t. to baseline of cell type **s** and the denominator corresponds to the sum of the proportions at baseline of the removed cell types_
$$
rv_{R,s} = {{\sum}_{i=1}^{1600} d_i\over {\sum}_{i=1}^{1600} {\sum}_{r∈R} B[i,r]}
$$