# This shell script executes R based deconvolution methods for removal of specific number of cell type removal scenario in parallel


# command line arguments provided to the scripts are,

# number of st datasets
# number of sc-ref datasets
# path to st dataset
# path to sc-ref dataset
# path to saving deconvolution results


z=$1
y=$2


for y in `seq 1 $y`
do
    for z in `seq $z $z`
    do
        nohup Rscript CARD.R $z $y $3 $4 $5 &
        nohup Rscript RCTD.R $z $y $3 $4 $5 &
        nohup Rscript Seurat.R $z $y $3 $4 $5 &
        nohup Rscript SPOTlight.R $z $y $3 $4 $5 &
        nohup Rscript SCDC.R $z $y $3 $4 $5 &
        nohup Rscript MuSiC.R $z $y $3 $4 $5 &
    done
done