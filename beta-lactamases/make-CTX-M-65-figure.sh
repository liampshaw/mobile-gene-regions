#!/bin/bash
#grep "CTX-M-65" ~/Downloads/test-bl/output/gene_variants/sequence_assignments/CTX-M-65.csv | cut -d ',' -f 1 > tmp_ctxm65.txt
#grep "CTX-M-27" ~/Downloads/test-bl/output/gene_variants/sequence_assignments/CTX-M-65.csv | cut -d ',' -f 1 > tmp_ctxm27.txt
#grep "CTX-M-14" ~/Downloads/test-bl/output/gene_variants/sequence_assignments/CTX-M-65.csv | cut -d ',' -f 1 > tmp_ctxm14.txt

grep "0 diffs"  ~/Downloads/test-bl/output/pangraph/CTX-M-65/CTX-M-65_extracted.fa | shuf | head -n 15 | cut -d ' ' -f 1 | tr -d '>' > tmp_ctxm65_0diff.txt
grep "1 diffs"  ~/Downloads/test-bl/output/pangraph/CTX-M-65/CTX-M-65_extracted.fa | cut -d ' ' -f 1 | tr -d '>' > tmp_ctxm65_1diff.txt
grep "CTX-M-24" ~/Downloads/test-bl/output/gene_variants/sequence_assignments/CTX-M-65.csv | cut -d ',' -f 1 > tmp_ctxm24.txt
grep -f tmp_ctxm65_1diff.txt tmp_ctxm24.txt > tmp_ctxm65_1diff_ctxm24.txt # there are 15 matches


grep "CTX-M-14" ~/Downloads/test-bl/output/gene_variants/sequence_assignments/CTX-M-65.csv | cut -d ',' -f 1 > tmp_ctxm14.txt
grep "2 diffs" ~/Downloads/test-bl/output/pangraph/CTX-M-65/CTX-M-65_extracted.fa | cut -d ' ' -f 1 |  tr -d '>'   > tmp_ctxm65_2diff.txt

# find those CTX-M-14 that have only 2 diffs, random subsample
grep -f tmp_ctxm65_2diff.txt tmp_ctxm14.txt | shuf | head -n 15 > tmp_ctxm65_2diff_ctxm14.txt

# plot linear blocks
python ../scripts/plot_linear_blocks_altair.py --block_csv ~/Downloads/test-bl/output/pangraph/CTX-M-65/CTX-M-65_pangraph.json.blocks.csv --output CTX-M-65-figure-v2/plot_ctxm65_0diff_subset.html --strain_list tmp_ctxm65_0diff.txt
python ../scripts/plot_linear_blocks_altair.py --block_csv ~/Downloads/test-bl/output/pangraph/CTX-M-65/CTX-M-65_pangraph.json.blocks.csv --output CTX-M-65-figure-v2/plot_ctxm65_1diff_ctxm24_subset.html --strain_list tmp_ctxm65_1diff_ctxm24.txt
python ../scripts/plot_linear_blocks_altair.py --block_csv ~/Downloads/test-bl/output/pangraph/CTX-M-65/CTX-M-65_pangraph.json.blocks.csv --output CTX-M-65-figure-v2/plot_ctxm65_2diff_ctxm14_subset.html --strain_list tmp_ctxm65_2diff_ctxm14.txt

