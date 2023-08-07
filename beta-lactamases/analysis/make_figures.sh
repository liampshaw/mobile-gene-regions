#!/bin/bash

./01_make_breakpoint_plots.sh  
./02_make_CTX-M-65_figure.sh  
Rscript 03_make_breakpoint_AUC_plots.R
./04_make_entropy_transposase_curves.sh
