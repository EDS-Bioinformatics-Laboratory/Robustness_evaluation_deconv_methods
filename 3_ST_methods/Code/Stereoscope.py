################################################################################
# Stereoscope posits a probabilistic model of spatial transcriptomics and an
# associated method for the deconvolution of cell type profiles using a
# single-cell RNA sequencing reference dataset.
#
# References
# Publication: https://www.nature.com/articles/s42003-020-01247-y
# Tutorial: https://docs.scvi-tools.org/en/stable/user_guide/models/stereoscope.html
# GitHub: https://github.com/almaan/stereoscope
#
# Implemented by Utkarsh Mahamune
# Bioinformatics Laboratory | Amsterdam UMC (Location AMC)
################################################################################


from matplotlib.figure import Figure
import numpy as np
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from scvi.external import RNAStereoscope, SpatialStereoscope

import warnings

warnings.filterwarnings("ignore")

import os
import sys

import jax

# scvi.settings.seed = 0

z = int(sys.argv[1])  # index of simulated ST dataset
y = int(sys.argv[2])  # number of single-cell reference datasets
arg3 = str(sys.argv[3])  # path to st data
arg4 = str(sys.argv[4])  # path to sc-ref data
arg5 = str(sys.argv[5])  # path to saving results


import time

startTime = time.time()

Data = "../../Data/Processed/"
Results = "../Results/"


for y in range(1, y + 1):
    sc_ref_data = sc.read_h5ad(arg4 + "sc.ref.data." + str(y) + ".h5ad")
    print(sc_ref_data)
    celltype_key = "blue.main"

    sc_ref_data.layers["counts"] = sc_ref_data.X.copy()
    # sc.pp.normalize_total(sc_ref_data, target_sum = 1e5)
    # sc.pp.log1p(sc_ref_data)
    # 
    # sc.pp.highly_variable_genes(sc_ref_data,
    #                             n_top_genes = 5000,
    #                             subset = True,
    #                             span = 1)

    z = int(sys.argv[1])
    for z in range(z, z + 1):
        spatial_data = sc.read_h5ad(arg3 + "spatial_data." + str(z) + ".h5ad")
        spatial_data.var_names_make_unique()

        # shared genes
        # intersect = np.intersect1d(list(sc_ref_data.var_names), list(spatial_data.var['features']))
        # 
        # sc_ref_data = sc_ref_data[:, intersect].copy()
        # spatial_data = spatial_data[:, intersect].copy()

        # Setup the scRNA-seq model
        
        RNAStereoscope.setup_anndata(sc_ref_data, layer = "counts",
                                      labels_key = celltype_key)

        # Train the scRNA-seq model
        sc_model = RNAStereoscope(sc_ref_data)
        sc_model.train(max_epochs=100)
        # sc_model.history["elbo_train"][10:].plot()
        # sc_model.save("stereo_sc_epochs.pdf", overwrite=True)
        
        # Infer proportions for spatial data
        spatial_data.layers["counts"] = spatial_data.X.copy()
        SpatialStereoscope.setup_anndata(spatial_data, layer="counts")
        
        # Train spatial data
        spatial_model = SpatialStereoscope.from_rna_model(spatial_data, sc_model)
        spatial_model.train(max_epochs=2000)
        # spatial_model.history["elbo_train"][10:].plot()
        # spatial_model.save("cell2loc_spatial_epochs.pdf", overwrite=True)

        # Stereoscope results
        spatial_data.obsm["deconvolution"] = spatial_model.get_proportions()
        
        print(spatial_model.get_proportions())
        
        spatial_model.get_proportions().to_csv(
            arg5
            + "Stereoscope/Stereoscope_result."
            + str(z)
            + "-"
            + str(y)
            + ".csv")


print("Execution time for stereoscope is: ", (time.time() - startTime) / 60)


