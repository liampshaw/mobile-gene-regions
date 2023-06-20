while read gene;
do
	echo $gene
	grep ">" /Users/Liam/Downloads/test-bl/output/pangraph/$gene/"$gene"_extracted.fa  | tr -d '>'  | cut -d ' ' -f 1 > breakpoint_plots/"$gene"_accs.txt
	python ../scripts/remove_duplicate_strains.py --strain_list breakpoint_plots/"$gene"_accs.txt --metadata ../data/metadata.csv --metadata_fields TaxGenus,Year,Country --output breakpoint_plots/"$gene"_accs_deduplicated.txt
	python ../scripts/plot_output_dists_altair.py --dist_csv /Users/Liam/Downloads/test-bl/output/pangraph/"$gene"/"$gene".output_dists.csv --gene_of_interest $gene --strain_list breakpoint_plots/"$gene"_accs_deduplicated.txt --output_html breakpoint_plots/"$gene"_deduplicated.html
        python ../scripts/plot_output_dists_altair.py --dist_csv /Users/Liam/Downloads/test-bl/output/pangraph/"$gene"/"$gene".output_dists.csv --gene_of_interest $gene --strain_list breakpoint_plots/"$gene"_accs.txt --output_html breakpoint_plots/"$gene"_all.html
	python ../scripts/calculate_AUC.py --dist_csv /Users/Liam/Downloads/test-bl/output/pangraph/"$gene"/"$gene".output_dists.csv --gene_of_interest $gene > AUCs/"$gene"_AUCs_focal_gene.txt
        python ../scripts/calculate_AUC.py --dist_csv /Users/Liam/Downloads/test-bl/output/pangraph/"$gene"/"$gene".output_dists.csv --gene_of_interest $gene --strain_list breakpoint_plots/"$gene"_accs_deduplicated.txt > AUCs/"$gene"_AUCs_focal_gene_deduplicated.txt
done < genes.txt

