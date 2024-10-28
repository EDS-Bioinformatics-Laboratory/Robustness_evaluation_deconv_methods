@echo off
REM Batch script for executing R-based deconvolution methods for a removal scenario

REM Command line arguments:
REM %1 - index of st datasets
REM %2 - number of sc-ref datasets
REM %3 - removal scenario

set z=%1
set y=%2
set st=../../2_Simulating_ST_data/Results/Spatial_Data/ REM Path to ST dataset
set sc_ref=../../1_Generate_sc_ref_data/Results/%3/    REM Path to sc-ref dataset
set decon_results=../Results/%3/                      REM Path to saving deconvolution results

FOR /L %%i IN (1,1,%y%) DO (
    FOR /L %%j IN (%z%,1,%z%) DO (
        start /B Rscript CARD.R %%i %%j %st% %sc_ref% %decon_results% > logs\CARD_%3_%%j_%%i.txt
        start /B Rscript RCTD.R %%i %%j %st% %sc_ref% %decon_results% > logs\RCTD_%3_%%j_%%i.txt
        start /B Rscript Seurat.R %%i %%j %st% %sc_ref% %decon_results% > logs\Seurat_%3_%%j_%%i.txt
        start /B Rscript SPOTlight.R %%i %%j %st% %sc_ref% %decon_results% > logs\SPOTlight_%3_%%j_%%i.txt
        start /B Rscript SCDC.R %%i %%j %st% %sc_ref% %decon_results% > logs\SCDC_%3_%%j_%%i.txt
        start /B Rscript MuSiC.R %%i %%j %st% %sc_ref% %decon_results% > logs\MuSiC_%3_%%j_%%i.txt
    )
)
