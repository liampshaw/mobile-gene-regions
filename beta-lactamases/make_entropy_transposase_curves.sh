while read gene;
do
	echo $gene
        grep ">" /Users/Liam/Downloads/test-bl/output/pangraph/$gene/"$gene"_extracted.fa  | tr -d '>'  | cut -d ' ' -f 1 > curves/"$gene"_accs_all.txt
	grep "0 diffs" /Users/Liam/Downloads/test-bl/output/pangraph/$gene/"$gene"_extracted.fa  | tr -d '>'  | cut -d ' ' -f 1 > curves/"$gene"_accs_focal_gene.txt
	python ../scripts/remove_duplicate_strains.py --strain_list curves/"$gene"_accs_focal_gene.txt --metadata ../data/metadata.csv --metadata_fields TaxGenus,Year,Country --output curves/"$gene"_accs_focal_gene_deduplicated.txt
	python ../scripts/remove_duplicate_strains.py --strain_list curves/"$gene"_accs_all.txt --metadata ../data/metadata.csv --metadata_fields TaxGenus,Year,Country --output curves/"$gene"_accs_all_deduplicated.txt
	# generate entropies for the focal gene subset
	python ../scripts/positional_entropy.py --json /Users/Liam/Downloads/test-bl/output/pangraph/"$gene"/"$gene"_pangraph.json --normalise --consensus --subset curves/"$gene"_accs_focal_gene.txt > entropies/"$gene"_entropies_focal_gene.txt 
	python ../scripts/positional_entropy.py --json /Users/Liam/Downloads/test-bl/output/pangraph/"$gene"/"$gene"_pangraph.json --normalise --consensus --subset curves/"$gene"_accs_all.txt > entropies/"$gene"_entropies_all.txt 
        python ../scripts/positional_entropy.py --json /Users/Liam/Downloads/test-bl/output/pangraph/"$gene"/"$gene"_pangraph.json --normalise --consensus --subset curves/"$gene"_accs_focal_gene_deduplicated.txt >    entropies/"$gene"_entropies_focal_gene_deduplicated.txt
	python ../scripts/positional_entropy.py --json /Users/Liam/Downloads/test-bl/output/pangraph/"$gene"/"$gene"_pangraph.json --normalise --consensus --subset curves/"$gene"_accs_all_deduplicated.txt > entropies/"$gene"_entropies_all_deduplicated.txt 
	# plot the transposase density
	python plot_transposase_density.py --dist_csv /Users/Liam/Downloads/test-bl/output/pangraph/"$gene"/"$gene".output_dists.csv --gff_file /Users/Liam/Downloads/test-bl/input/gffs/"$gene"_annotations.gff --gene_name "$gene" --output_html curves/"$gene"-curves-focal-gene.html --flanking_width 5000 --entropies entropies/"$gene"_entropies_focal_gene.txt --strain_list curves/"$gene"_accs_focal_gene.txt
        python plot_transposase_density.py --dist_csv /Users/Liam/Downloads/test-bl/output/pangraph/"$gene"/"$gene".output_dists.csv --gff_file /Users/Liam/Downloads/test-bl/input/gffs/"$gene"_annotations.gff --gene_name "$gene" --output_html curves/"$gene"-curves-focal-genededuplicated.html --flanking_width  5000 --entropies entropies/"$gene"_entropies_focal_gene_deduplicated.txt --strain_list curves/"$gene"_accs_focal_gene_deduplicated.txt
        python plot_transposase_density.py --dist_csv /Users/Liam/Downloads/test-bl/output/pangraph/"$gene"/"$gene".output_dists.csv --gff_file /Users/Liam/Downloads/test-bl/input/gffs/"$gene"_annotations.gff --gene_name "$gene" --output_html curves/"$gene"-curves-all-deduplicated.html --flanking_width  5000 --entropies entropies/"$gene"_entropies_all_deduplicated.txt --strain_list curves/"$gene"_accs_all_deduplicated.txt
        python plot_transposase_density.py --dist_csv /Users/Liam/Downloads/test-bl/output/pangraph/"$gene"/"$gene".output_dists.csv --gff_file /Users/Liam/Downloads/test-bl/input/gffs/"$gene"_annotations.gff --gene_name "$gene" --output_html curves/"$gene"-curves-all.html --flanking_width  5000 --entropies entropies/"$gene"_entropies_all.txt --strain_list curves/"$gene"_accs_all.txt
done < genes.txt

