# See configfile in e.g. configs/laptop_config.yaml for parameters
DB = config["DB"]
# #FOCAL_GENE_DICT = config["focal_gene_dict"]
# FOCAL_GENE_DICT = {"CMY-2": "CMY", 
# 					"CTX-M-15": "CTX-M",
# 					"CTX-M-65":"CTX-M",
# 					"GES-24":"GES",
# 					"IMP-4":"IMP",
# 					"KPC-2": "KPC",
# 					"NDM-1":"NDM",
# 					"PER-1":"PER",  
# 					"VIM-1":"VIM",
# 					"TEM-1": "TEM",
# 					"OXA-10": "OXA",
# 					"OXA-48":"OXA"}

FOCAL_GENE_DICT = {"GES-24": "GES"}


rule run_pangraph:
	input:
		expand("output/pangraph/{focal_gene}/{focal_gene}_pangraph.json", focal_gene=FOCAL_GENE_DICT.keys())

rule calculate_distances:
	input:
		expand("output/pangraph/{focal_gene}/{focal_gene}.output_dists.csv", focal_gene=FOCAL_GENE_DICT.keys())

rule make_plots:
	input:
		expand("output/pangraph/{focal_gene}/{focal_gene}_pangraph.json", focal_gene=FOCAL_GENE_DICT.keys())




#############################
# GETTING CONTIG DATA READY #
#############################
rule combine_fastas_containing_gene:
	input:
		"output/analysis/{focal_gene}_accessions.txt"
	output:
		"output/contigs/{focal_gene}_combined_contigs.fa"
	run:
		shell("while read f; do cat "+f"{config['fastadir']}"+"/$f.fa >> {output}; done < {input}")

rule extract_genes_from_contigs:
	input:
		gene_fasta="data/focal_genes/{focal_gene}.fa",
		input_fasta= lambda wildcards:f"output/contigs/{FOCAL_GENE_DICT[wildcards.focal_gene]}"+"_combined_contigs.fa"
	params:
		snp_threshold=int(config["snp_threshold"])
	output:
		"output/analysis/sequence/{focal_gene}_seqs_extracted_from_contigs.fa"
	shell:
		"python scripts/extract_region_around_gene.py --gene {input.gene_fasta} \
												--input_fasta {input.input_fasta} \
												--output_fasta {output} \
												--upstream 0 \
												--downstream 0 \
												--threshold {params.snp_threshold}"

rule name_observed_genes:
	input:
		fasta="output/analysis/sequence/{focal_gene}_seqs_extracted_from_contigs.fa",
		variants="DB/variant_fasta/{focal_gene}.fa"
	output:
		"output/analysis/sequence_assignments/{focal_gene}.csv"
	shell:
		"python scripts/name_variants.py --variant_fasta {input.variants} --output_file {output} --input_fasta {input.fasta}"

#############################
# RUNNING PANGRAPH PIPELINE #
#############################

rule extract_region_around_focal_gene:
	input:
		input_fasta=lambda wildcards:f"output/contigs/{FOCAL_GENE_DICT[wildcards.focal_gene]}"+"_combined_contigs.fa"
	params:
		focal_gene="data/focal_genes/{focal_gene}.fa",
		prefix="output/pangraph/{focal_gene}/{focal_gene}_extracted",
		upstream=config["region_upstream"],
		downstream=config["region_downstream"],
		threshold=int(config["snp_threshold"])
	output:
		"output/pangraph/{focal_gene}/{focal_gene}_extracted.fa"#,
	shell: 
		"python scripts/extract_region_around_gene.py --gene {params.focal_gene} --input {input.input_fasta} \
		--upstream {params.upstream} --downstream {params.upstream} --complete --output_fasta {params.prefix}.fa\
		--threshold {params.threshold}"

rule calculate_snp_dists_extracted_seqs:
	input:
		"output/analysis/sequence/{focal_gene}_seqs_extracted_from_contigs.fa"
	output:
		"output/analysis/sequence/{focal_gene}_seqs_extracted_from_contigs.snps.tsv",
		"output/analysis/sequence/{focal_gene}_seqs_extracted_from_contigs.fa.dedup.aln",
		"output/analysis/sequence/{focal_gene}_seqs_extracted_from_contigs.fa.dedup.txt"
	run:
		shell("mafft --auto {input} > {input}.aln"),
		shell("snp-dists -q -m {input}.aln > {output}")
		shell("seqkit rmdup -s < {input}.aln > {input}.dedup.aln -D {input}.dedup.txt")

rule tree_for_focal_gene:
	input:
		"output/analysis/sequence/{focal_gene}_seqs_extracted_from_contigs.fa.dedup.aln"
	output:
		"output/analysis/sequence/{focal_gene}_seqs_extracted_from_contigs.fa.dedup.aln.refined.tre"
	run:
		shell("FastTree -quiet -nt -gtr {input} > {input}.tre")
		shell("python scripts/refine_tree.py --aln {input} --tree_in {input}.tre --tree_out {output}")


rule build_pangraph:
	input:
		"output/pangraph/{focal_gene}/{focal_gene}_extracted.fa"
	params:
		aligner=config["pangraph_aligner"],
		minlength=config["pangraph_minblocklength"]
	output:
		"output/pangraph/{focal_gene}/{focal_gene}_pangraph.json"
	shell:
		"pangraph build -k {params.aligner} --len {params.minlength} {input} > {output}" if config["pangraph_polish"]==False 
		else "pangraph build -k {params.aligner} --len {params.minlength} {input} | pangraph polish > {output}"


rule export_pangraph:
	input:
		"output/pangraph/{focal_gene}/{focal_gene}_pangraph.json"
	params:
		prefix="{focal_gene}_pangraph",
		outdir="output/pangraph/{focal_gene}/",
		edgeminlength=config["pangraph_edgeminlength"]
	output:
		"output/pangraph/{focal_gene}/{focal_gene}_pangraph.gfa",
		"output/pangraph/{focal_gene}/{focal_gene}_pangraph.fa"
	shell:
		"pangraph export --edge-minimum-length {params.edgeminlength} {input} \
							-p {params.prefix} \
							-o {params.outdir}" if config["panx_export"]==False else
		"pangraph export --edge-minimum-length {params.edgeminlength} {input} \
							-p {params.prefix} \
							-o {params.outdir} \
							--export-panX"


rule convert_pangraph_to_block_list:
	input:
		json="output/pangraph/{focal_gene}/{focal_gene}_pangraph.json",
		gfa="output/pangraph/{focal_gene}/{focal_gene}_pangraph.gfa"
	output:
		"output/pangraph/{focal_gene}/{focal_gene}_pangraph.gfa.coloured.gfa",
		"output/pangraph/{focal_gene}/{focal_gene}_pangraph.json.blocks.csv"
	shell:
		"python scripts/convert_pangraph_to_block_list.py --json {input.json} --gfa {input.gfa}"

rule find_focal_gene_block:
	input:
		focal_gene="data/focal_genes/{focal_gene}.fa",
		pangraph_fasta="output/pangraph/{focal_gene}/{focal_gene}_pangraph.fa"
	output:
		"output/pangraph/{focal_gene}/{focal_gene}.gene_block.txt"
	run:
		shell("makeblastdb -in {input.pangraph_fasta} -dbtype 'nucl'"),
		shell("blastn -query {input.focal_gene} -db {input.pangraph_fasta} -outfmt 6 | cut -f 2 > {output}")

rule compute_distances:
	input:
		block_csv="output/pangraph/{focal_gene}/{focal_gene}_pangraph.json.blocks.csv",
		gene_block="output/pangraph/{focal_gene}/{focal_gene}.gene_block.txt",
		snps="output/analysis/sequence/{focal_gene}_seqs_extracted_from_contigs.snps.tsv"
	output:
		"output/pangraph/{focal_gene}/{focal_gene}.output_dists.csv"
	shell:
		"python scripts/compute_distances.py --block_csv {input.block_csv} \
											--gene_block_file {input.gene_block} \
											--snps {input.snps}\
											--output {output}"
