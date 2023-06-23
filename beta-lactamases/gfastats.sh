while read gene;
do
	/Users/Liam/github/gfastats/build/bin/gfastats /Users/Liam/Downloads/test-bl/output/pangraph/"$gene"/"$gene"_pangraph.gfa > gfastats/"$gene"_gfastats.txt
	seqkit stats /Users/Liam/Downloads/test-bl/output/pangraph/"$gene"/"$gene"_pangraph.fa > gfastats/"$gene"_pangraph_size.txt
done < genes.txt
