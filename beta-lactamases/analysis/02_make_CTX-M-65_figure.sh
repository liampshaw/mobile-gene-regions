#!/bin/bash

mkdir -p CTX-M-65-figure

grep "0 diffs"  ../output/CTX-M-65/pangraph/CTX-M-65_extracted.fa | shuf | head -n 15 | cut -d ' ' -f 1 | tr -d '>' > tmp_ctxm65_0diff.txt
grep "1 diffs"  ../output/CTX-M-65/pangraph/CTX-M-65_extracted.fa | cut -d ' ' -f 1 | tr -d '>' > tmp_ctxm65_1diff.txt
grep "CTX-M-24" ../output/CTX-M-65/gene_diversity/CTX-M-65.csv | cut -d ',' -f 1 > tmp_ctxm24.txt
grep -f tmp_ctxm65_1diff.txt tmp_ctxm24.txt > tmp_ctxm65_1diff_ctxm24.txt # there are 15 matches


grep "CTX-M-14" ../output/CTX-M-65/gene_diversity/CTX-M-65.csv | cut -d ',' -f 1 > tmp_ctxm14.txt
grep "2 diffs" ../output/CTX-M-65/pangraph/CTX-M-65_extracted.fa | cut -d ' ' -f 1 |  tr -d '>'   > tmp_ctxm65_2diff.txt

# find those CTX-M-14 that have only 2 diffs, random subsample
grep -f tmp_ctxm65_2diff.txt tmp_ctxm14.txt | shuf | head -n 15 > tmp_ctxm65_2diff_ctxm14.txt

# plot linear blocks
python ../../scripts/plot_linear_blocks_altair.py --block_csv ../output/CTX-M-65/pangraph/CTX-M-65_pangraph.json.blocks.csv --output CTX-M-65-figure/plot_ctxm65_0diff_subset.html --strain_list tmp_ctxm65_0diff.txt
python ../../scripts/plot_linear_blocks_altair.py --block_csv ../output/CTX-M-65/pangraph/CTX-M-65_pangraph.json.blocks.csv --output CTX-M-65-figure/plot_ctxm65_1diff_ctxm24_subset.html --strain_list tmp_ctxm65_1diff_ctxm24.txt
python ../../scripts/plot_linear_blocks_altair.py --block_csv ../output/CTX-M-65/pangraph/CTX-M-65_pangraph.json.blocks.csv --output CTX-M-65-figure/plot_ctxm65_2diff_ctxm14_subset.html --strain_list tmp_ctxm65_2diff_ctxm14.txt

# Deduplicated version
python ../../scripts/plot_output_dists_altair.py --dist_csv ../output/CTX-M-65/pangraph/CTX-M-65.output_dists.csv --output_html CTX-M-65-figure/figS1_deduplicated.html --gene_of_interest CTX-M-65 --strain_list breakpoint_plots/CTX-M-65_accs_deduplicated.txt 
