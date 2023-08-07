#!/bin/bash

echo "############################"
echo "Setting up breakpoint plots..."
echo "############################"
./01_make_breakpoint_plots.sh  



echo "############################"
echo "Making CTX-M-65 figure..."
echo "############################"
./02_make_CTX-M-65_figure.sh  


echo "############################"
echo "Making breakpoint AUC plots..."
echo "############################"
Rscript 03_make_breakpoint_AUC_plots.R

echo "############################"
echo "Making annotation/entropy/breakpoint curves..."
echo "############################"
./04_make_entropy_transposase_curves.sh
