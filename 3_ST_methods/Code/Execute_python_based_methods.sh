# This script executes python based deconvolution methods for removal of specific number of cell type/s from reference dataset

# command line arguments provided to the scripts are,

# index of st datasets
# number of sc-ref datasets
# removal scenario

z=$1
y=$2
st="../../2_Simulating_ST_data/Results/Spatial_Data/" # path to st dataset
sc_ref="../../1_Generate_sc_ref_data/Results/$3/" # path sc reference dataset
decon_results="../Results/$3/" # path to saving deconvolution results

x="${3}_${y}_${z}"
nohup python3 Cell2Location.py $1 $2 $st $sc_ref $decon_results > logs/Cell2location_$x.txt &
nohup python3 Stereoscope.py $1 $2 $st $sc_ref $decon_results > logs/Stereoscope_$x.txt &
