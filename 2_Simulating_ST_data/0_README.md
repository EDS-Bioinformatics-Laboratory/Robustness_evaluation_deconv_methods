# Generating simulated spatial transcriptomics datasets


**Operating System(s) / version(s) used during development (and testing):**

MacOS Ventura 13.4.1


**Software environment:** 

- R version : 4.1.2 (2021-11-01)
- Rstudio IDE: RStudio 2023.03.0+386 "Cherry Blossom"



####**Conceptual description of methodology**


We generate a single spot and apply the same algorithm to the rest of the spots. In our simulated spatial transcriptomics, we consider each spot an independent entity, inferring no biological relevance among neighboring spots. This assumption suits our goal of establishing a ground truth and only focuses on the deconvolution of a spot for cell type proportion using the reference scRNA-seq data.

*  First we draw the number of cell types *C<sub>i</sub>* in spot *i* from a discrete uniform distribution *U* {*C*<sub>min</sub>,*C*<sub>max</sub>} and then randomly sample *C<sub>i</sub>* cell types.

*  We then determine the proportion of each cell type present in a spot by *C<sub>i</sub>* independent draws from a continuous uniform distribution *U*(0,1). The resulting values are normalized such that their sum equals 1. The resulting proportions *P<sub>i,k</sub>* of cell type *k* serve as the **ground truth** for our benchmarking study.

* Subsequently, the number of cells *D<sub>i</sub>* in a spot is drawn from a discrete uniform distribution *U* {*D*<sub>min</sub>,*D*<sub>max</sub>}.

* Next, using both the cell type proportions *P<sub>i,k</sub>* and the number of cells *D<sub>i </sub>*, we determine the number of cells per cell type *W<sub>i,k</sub>* by sampling from a multinomial distribution Multinomial(*D<sub>i</sub>* ; *P<sub>i,1</sub>,…, P<sub>i,C</sub>*).

* We then simulate the mRNA counts *M<sub>i,g</sub>* for gene *g* using the following equation:		
$$
{M_{i, g}} = \sum_{k = 1}^{C}{{R_{g, k}} * W_{i, k}}
$$
Here, *R<sub>g,k</sub>* denotes the mRNA count for gene *g* in cell type *k*. *R<sub>g,k</sub>* is determined per spot by calculating the ceiling of the average of the counts for gene *g* in five cells of cell type *k* randomly drawn from the reference dataset. Finally, to make the total number of counts per spot similar to what is observed in 10x Visium spatial transcriptomics data for those spots with more than 30,000 counts, we draw a random number from a discrete uniform distribution U{25000, 30000} and scale the per gene counts accordingly (function ‘downsampleMatrix’, package scuttle).

* We also generate dummy spatial coordinates for each spot in 40 by 40 layout for those deconvolution methods that require the spot coordinates as input. 


