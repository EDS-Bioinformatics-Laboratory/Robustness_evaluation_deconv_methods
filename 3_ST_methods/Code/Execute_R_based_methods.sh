# This shell script executes R-based deconvolution methods for a removal scenario in parallel for all the reference datasets provided through command line


# command line arguments provided to the scripts are,

# index of st datasets
# number of sc-ref datasets
# reomoval scenario


z=$1
y=$2
st="../../2_Simulating_ST_data/Results/Spatial_Data/" # path to st dataset
sc_ref="../../1_Generate_sc_ref_data/Results/$3/" # path sc reference dataset
decon_results="../Results/$3/" # path to saving deconvolution results



for y in `seq 1 $y`
do
    for z in `seq $z $z`
    do
        x="${3}_${y}_${z}"
        nohup Rscript CARD.R $z $y $st $sc_ref $decon_results > logs/CARD_$x.txt &
        nohup Rscript RCTD.R $z $y $st $sc_ref $decon_results > logs/RCTD_$x.txt &
        nohup Rscript Seurat.R $z $y $st $sc_ref $decon_results > logs/Seurat_$x.txt&
        nohup Rscript SPOTlight.R $z $y $st $sc_ref $decon_results > logs/SPOTlight_$x.txt &
        nohup Rscript SCDC.R $z $y $st $sc_ref $decon_results > logs/SCDC_$x.txt &
        nohup Rscript MuSiC.R $z $y $st $sc_ref $decon_results > logs/MuSiC_$x.txt &
    done
done