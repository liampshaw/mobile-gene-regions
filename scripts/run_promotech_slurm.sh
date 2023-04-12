#!/bin/bash

# NOTE: not part of the pipeline but required to generate locations of 
# 40nt stretches of contigs in database that contain a promoter with a 
# probability > 0.6 as determined by PromoTech
# Date: 29/03/2023

module load Anaconda3/2020.11
eval "$(conda shell.bash hook)"
conda activate promotech_env
overall_outdir=promoter_predictions
outdir=tmp-large-$SLURM_ARRAY_TASK_ID
fasta=$(sed "${SLURM_ARRAY_TASK_ID}q;d" fastas.txt)
name=$(echo $fasta | rev | cut -d '.' -f 2- | cut -d '/' -f 1 | rev )
echo $name
time python /well/shaw/users/amu125/programs/PromoTech/promotech.py -pg -m RF-HOT -f $fasta -o $outdir
echo "python /well/shaw/users/amu125/programs/PromoTech/promotech.py -pg -m RF-HOT -f $fasta -o $outdir"

time python /well/shaw/users/amu125/programs/PromoTech/promotech.py -g -t 0.6 -i $outdir -o $outdir
echo "python /well/shaw/users/amu125/programs/PromoTech/promotech.py -g -t 0.6 -i $outdir -o $outdir"
rm $outdir/*data
cp $outdir/genome_predictions.csv $overall_outdir/$name-promoter_predictions.tsv
rm -r $outdir

