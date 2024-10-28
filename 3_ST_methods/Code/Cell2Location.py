################################################################################
# Cell2location is a principled Bayesian model that can resolve fine-grained
# cell types in spatial transcriptomic data and create comprehensive cellular
# maps of diverse tissues.
#
# References
# Publication: https://www.nature.com/articles/s41587-021-01139-4
# Tutorial: https://cell2location.readthedocs.io/en/latest/
# GitHub: https://github.com/BayraktarLab/cell2location
#
# Implemented by Utkarsh Mahamune
# Bioinformatics Laboratory | Amsterdam UMC (Location AMC)
################################################################################

import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import subprocess
import os
import cell2location
from matplotlib import rcParams

rcParams["pdf.fonttype"] = 42  # enables correct plotting of text
import seaborn as sns
from scipy.sparse import csr_matrix
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel

import sys
import jax

# scvi.settings.seed = 0

z = int(sys.argv[1])  # index of simulated ST dataset
y = int(sys.argv[2])  # number of single-cell reference datasets
arg3 = str(sys.argv[3])  # path to st data
arg4 = str(sys.argv[4])  # path to sc-ref data
arg5 = str(sys.argv[5])  # path to saving results

# value for hyper parameter
if (z == 3):
  num_cells = 5
else:
    num_cells = 13

import time

startTime = time.time()

Data = "../../Data/Processed/"
Results = "../Results/"

for y in range(1, y + 1):
    # loading in the sc ref and spatial data
    sc_ref_data = sc.read_h5ad(arg4 + "sc.ref.data." + str(y) + ".h5ad")
    celltype_key = "blue.main"

    sc_ref_data.X = csr_matrix(sc_ref_data.X)

    # # scanpy filtering
    # sc.pp.filter_genes(sc_ref_data, min_cells = 0)
    # sc.pp.filter_cells(sc_ref_data, min_genes = 0)

    # # cell2location filtering of genes
    # selected = filter_genes(
    #     sc_ref_data,
    #     cell_count_cutoff = 5,
    #     cell_percentage_cutoff2 = 0.03,
    #     nonz_mean_cutoff = 2,
    # )
    # sc_ref_data = sc_ref_data[:, selected].copy()

    # preparing anndata for sc ref data
    cell2location.models.RegressionModel.setup_anndata(
        sc_ref_data, labels_key=celltype_key
    )

    # create and train model for sc ref data
    mod = RegressionModel(sc_ref_data)
    mod.train(max_epochs = 1000, batch_size = 2500, train_size = 1, lr = 0.002, use_gpu = True)
    # mod.plot_history(10)
    # plt.savefig("cell2loc_sc_epochs.pdf")
        
    # export estimated cell abundance (summary of posterior distribution
    sc_ref_data = mod.export_posterior(
        sc_ref_data,
        sample_kwargs={"num_samples": 1000, "batch_size": 2500, "use_gpu": True},
    )

    # extracting reference cell types signatures as a dataframe from binomial regression model
    if "means_per_cluster_mu_fg" in sc_ref_data.varm.keys():
        inf_aver = sc_ref_data.varm["means_per_cluster_mu_fg"][
            [
                f"means_per_cluster_mu_fg_{i}"
                for i in sc_ref_data.uns["mod"]["factor_names"]
            ]
        ].copy()
    else:
        inf_aver = sc_ref_data[
            [
                f"means_per_cluster_mu_fg_{i}"
                for i in sc_ref_data.uns["mod"]["factor_names"]
            ]
        ].copy()
    inf_aver.columns = sc_ref_data.uns["mod"]["factor_names"]

    z = int(sys.argv[1])
    for z in range(z, z + 1):
        spatial_data = sc.read_h5ad(arg3 + "spatial_data." + str(z) + ".h5ad")
        spatial_data.var_names = spatial_data.var['features']

        spatial_data.X = csr_matrix(spatial_data.X)

        # finding shared genes
        intersect = np.intersect1d(spatial_data.var_names, inf_aver.index)
        spatial_data = spatial_data[:, intersect].copy()
        inf_aver = inf_aver.loc[intersect, :].copy()

        # preparing anndata for spatial data
        cell2location.models.Cell2location.setup_anndata(spatial_data)

        # create and train model for spatial data
        mod = cell2location.models.Cell2location(
            spatial_data,
            cell_state_df = inf_aver,
            # the expected average cell abundance: tissue-dependent
            # hyper-prior which can be estimated from paired histology:
            N_cells_per_location = num_cells,
            # hyperparameter controlling normalisation of
            # within-experiment variation in RNA detection:
            # suggested by developer
            detection_alpha = 200
        )
        mod.train(max_epochs = 10000, batch_size = None, train_size = 1, use_gpu = True)
        # mod.plot_history(1000)
        # plt.savefig("cell2loc_spatial_epochs.pdf")

        # export estimated cell abundance (summary of posterior distribution)
        spatial_data = mod.export_posterior(
            spatial_data,
            sample_kwargs={
                "num_samples": 1000,
                "batch_size": mod.adata.n_obs,
                "use_gpu": True,
            },
        )

        result_cell2location = pd.DataFrame(
            spatial_data.obsm["q05_cell_abundance_w_sf"]
        )
        result_cell2location.to_csv(
            arg5 + "Cell2Location/Cell2Location." + str(z) + "-" + str(y) + ".csv"
        )

print("Execution time for cell2location is: ", (time.time() - startTime) / 60)


