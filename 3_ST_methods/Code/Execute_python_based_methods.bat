@echo off
REM This script executes Python-based deconvolution methods for removal of specific number of cell types from the reference dataset

REM Command line arguments provided to the scripts are:
REM 1. Index of ST datasets
REM 2. Number of SC-ref datasets
REM 3. Removal scenario

SET "st=..\..\2_Simulating_ST_data\Results\Spatial_Data\"  REM Path to ST dataset
SET "sc_ref=..\..\1_Generate_sc_ref_data\Results\%3\"      REM Path to SC reference dataset
SET "decon_results=..\Results\%3\"                         REM Path to saving deconvolution results

START /B python Cell2Location.py %1 %2 %st% %sc_ref% %decon_results% > logs\Cell2location_%3_%1_%2.txt
START /B python Stereoscope.py %1 %2 %st% %sc_ref% %decon_results% > logs\Stereoscope_%3_%1_%2.txt
