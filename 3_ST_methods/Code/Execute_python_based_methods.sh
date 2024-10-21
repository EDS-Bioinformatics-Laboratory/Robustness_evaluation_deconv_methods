# This script executes python based deconvolution methods for removal of specific number of cell type/s from reference dataset

# command line arguments provided to the scripts are,

# number of st datasets
# number of sc-ref datasets
# path to st data
# path to sc-ref data
# path to saving results

nohup python3 Cell2Location.py $1 $2 $3 $4 $5 &
nohup python3 Stereoscope.py $1 $2 $3 $4 $5 &