#!/bin/bash

grep "# contigs" *txt | cut -d ':' -f 3 | tr -d ' ' > stats_contigs
grep "# bubbles" *txt | cut -d ':' -f 3 | tr -d ' ' > stats_bubbles
grep "Contig N50" *txt | cut -d ':' -f 3 | tr -d ' ' > stats_contig_N50
grep "circular segments" *txt | cut -d ':' -f 3 | tr -d ' ' > stats_circular_segments
grep "# scaffolds" *txt | cut -d ':' -f 3 | tr -d ' ' > stats_scaffolds
grep "Total contig length" *txt | cut -d ':' -f 3 | tr -d ' ' > stats_total_contig_length
grep "Average contig length" *txt | cut -d ':' -f 3 | tr -d ' ' > stats_average_contig_length
tail -n 1 *_pangraph_size.txt | awk '{print $5}' | tr -d ',' | awk 'NF' > stats_fasta

echo "gene,sequences,blocks,total_size,average_block_length,block_N50,bubbles,circular_segments,total_fasta_size" > combined_stats.csv
paste -d ',' ../genes.txt stats_scaffolds stats_contigs stats_total_contig_length stats_average_contig_length stats_contig_N50 stats_bubbles stats_circular_segments stats_fasta >> combined_stats.csv

