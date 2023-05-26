# See configfile in e.g. configs/laptop_config.yaml for parameters

FOCAL_GENES = ["mcr-1.1"]
# you have two input files:
# focal gene: input/focal_genes/mcr-1.1.fa
# contigs: input/contigs/mcr1.1_contigs.fa 
# (optional gff: input/gffs/mcr1.1_annotations.gff 
# if using include_gff: True in config file)

# N.B. if trying to get named variants from CARD DB this may break 
# if gene name is not in CARD / does not match 

DB = config["DB"]

rule prepare_DB:
	input:
		expand("output/gene_variants/{gene}.fa", gene=FOCAL_GENES), # Removed for now
		expand("output/gene_variants/sequence_assignments/{gene}.csv", gene=FOCAL_GENES)

rule run_pangraph:
	input:
		expand("output/pangraph/{gene}/{gene}_pangraph.json", gene=FOCAL_GENES)

rule calculate_distances:
	input:
		expand("output/pangraph/{gene}/{gene}.output_dists.csv", gene=FOCAL_GENES)

rule make_plots:
	input:
		expand("output/pangraph/{gene}/plots/linear_blocks.pdf", gene=FOCAL_GENES),
		expand("output/pangraph/{gene}/plots/{gene}_breakpoint_distances-all.pdf", 
			gene=FOCAL_GENES),
		expand("output/pangraph/{gene}/plots/bandage.log_file", gene=FOCAL_GENES),
		expand("output/pangraph/{gene}/plots/{gene}_positional_entropies_consensus_relative.pdf", gene=FOCAL_GENES),
		expand("output/pangraph/{gene}/plots/NJ_tree_central_gene.pdf", gene=FOCAL_GENES),
		expand("output/pangraph/{gene}/plots/{gene}_linear_blocks.html", gene=FOCAL_GENES),
		expand("output/pangraph/{gene}/plots/{gene}_linear_blocks_deduplicated.html", gene=FOCAL_GENES),
		expand("output/pangraph/{gene}/plots/{gene}_ecdf.html", gene=FOCAL_GENES),
		expand("output/pangraph/{gene}/plots/bandage.log_file", gene=FOCAL_GENES)

# # Unsure whether to include this
rule extract_genes_DB:
	input:
		"input/focal_genes/{gene}.fa"
	params:
		DB="data/CARD_db.fa"
	output:
		"output/gene_variants/{gene}.fa"
	shell:
		"python scripts/extract_gene_DB.py {params.DB} {input} CARD {output}"

rule extract_genes_from_contigs:
	input:
		gene_fasta="input/focal_genes/{gene}.fa",
		input_fasta="input/contigs/{gene}_contigs.fa"
	params:
		snp_threshold=int(config["snp_threshold"])
	output:
		"output/gene_variants/sequence/{gene}_seqs.fa"
	shell:
		"python scripts/extract_region_around_gene.py --gene {input.gene_fasta} \
												--input_fasta {input.input_fasta} \
												--output_fasta {output} \
												--upstream 0 \
												--downstream 0 \
												--threshold {params.snp_threshold}"

rule assign_variants:
	input:
		fasta="output/gene_variants/sequence/{gene}_seqs.fa",
		variants="output/gene_variants/{gene}.fa"
	output:
		"output/gene_variants/sequence_assignments/{gene}.csv"
	shell:
		"python scripts/name_variants.py --variant_fasta {input.variants} --output_file {output} --input_fasta {input.fasta}"

#############################
# RUNNING PANGRAPH PIPELINE #
#############################

rule extract_region_around_focal_gene:
	input:
		input_fasta="input/contigs/{gene}_contigs.fa"
	params:
		gene="input/focal_genes/{gene}.fa",
		prefix="output/pangraph/{gene}/{gene}_extracted",
		upstream=config["region_upstream"],
		downstream=config["region_downstream"],
		threshold=int(config["snp_threshold"])
	output:
		"output/pangraph/{gene}/{gene}_extracted.fa"
	shell: 
		"python scripts/extract_region_around_gene.py --gene {params.gene} --input {input.input_fasta} \
		--upstream {params.upstream} --downstream {params.upstream} --complete --output_fasta {params.prefix}.fa\
		--threshold {params.threshold}"

rule calculate_snp_dists_extracted_seqs:
	input:
		"output/gene_variants/sequence/{gene}_seqs.fa"
	output:
		"output/gene_variants/sequence/{gene}_seqs.snps.tsv",
		"output/gene_variants/sequence/{gene}_seqs.fa.dedup.aln",
		"output/gene_variants/sequence/{gene}_seqs.fa.dedup.txt"
	run:
		shell("mafft --auto {input} > {input}.aln"),
		shell("snp-dists -q -m {input}.aln > {output}")
		shell("seqkit rmdup -s < {input}.aln > {input}.dedup.aln -D {input}.dedup.txt")

rule tree_for_focal_gene:
	input:
		"output/gene_variants/sequence/{gene}_seqs.fa.dedup.aln"
	output:
		"output/gene_variants/sequence/{gene}_seqs.fa.dedup.aln.refined.tre"
	run:
		shell("FastTree -quiet -nt -gtr {input} > {input}.tre")
		shell("sed -e 's/:.*//g' {input} > {input}.renamed")
		shell("python scripts/refine_tree.py --aln {input}.renamed --tree_in {input}.tre --tree_out {output}")


rule build_pangraph:
	input:
		"output/pangraph/{gene}/{gene}_extracted.fa"
	params:
		aligner=config["pangraph_aligner"],
		minlength=config["pangraph_minblocklength"]
	output:
		"output/pangraph/{gene}/{gene}_pangraph.json"
	shell:
		"pangraph build -k {params.aligner} --len {params.minlength} {input} > {output}" if config["pangraph_polish"]==False 
		else "pangraph build -k {params.aligner} --len {params.minlength} {input} | pangraph polish > {output}"


rule export_pangraph:
	input:
		"output/pangraph/{gene}/{gene}_pangraph.json"
	params:
		prefix="{gene}_pangraph",
		outdir="output/pangraph/{gene}/",
		edgeminlength=config["pangraph_edgeminlength"]
	output:
		"output/pangraph/{gene}/{gene}_pangraph.gfa",
		"output/pangraph/{gene}/{gene}_pangraph.fa"
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
		json="output/pangraph/{gene}/{gene}_pangraph.json",
		gfa="output/pangraph/{gene}/{gene}_pangraph.gfa"
	output:
		"output/pangraph/{gene}/{gene}_pangraph.gfa.coloured.gfa",
		"output/pangraph/{gene}/{gene}_pangraph.json.blocks.csv"
	shell:
		"python scripts/convert_pangraph_to_block_list.py --json {input.json} --gfa {input.gfa}"

rule find_focal_gene_block:
	input:
		gene="input/focal_genes/{gene}.fa",
		pangraph_fasta="output/pangraph/{gene}/{gene}_pangraph.fa"
	output:
		"output/pangraph/{gene}/{gene}.gene_block.txt"
	run:
		shell("makeblastdb -in {input.pangraph_fasta} -dbtype 'nucl'"),
		shell("blastn -query {input.gene} -db {input.pangraph_fasta} -outfmt 6 | cut -f 2 > {output}")


rule gene_locations:
	input:
		gene="input/focal_genes/{gene}.fa",
		db="output/pangraph/{gene}/{gene}_extracted.fa",
		gene_block="output/pangraph/{gene}/{gene}.gene_block.txt"
	output:
		"output/pangraph/{gene}/{gene}_gene_locations_block.txt"
	run:
		shell("makeblastdb -in {input.db} -dbtype 'nucl'")
		shell("blastn -max_hsps 10000 -query {input.gene} -db {input.db} -outfmt '6 sseqid sstart send' > {output}")

rule compute_distances:
	input:
		gene_fasta="input/focal_genes/{gene}.fa",
		block_csv="output/pangraph/{gene}/{gene}_pangraph.json.blocks.csv",
		gene_block="output/pangraph/{gene}/{gene}.gene_block.txt",
		snps="output/gene_variants/sequence/{gene}_seqs.snps.tsv",
		#locations="output/pangraph/{gene}/{gene}_gene_locations.txt",
		pangraph_fasta="output/pangraph/{gene}/{gene}_pangraph.fa"
	output:
		gene_block_fasta="output/pangraph/{gene}/{gene}.gene_block.fa",
		gene_rel_locations="output/pangraph/{gene}/{gene}_relative_locations_block.txt",
		dists="output/pangraph/{gene}/{gene}.output_dists.csv"
	run:
		shell("./scripts/extract_seq_from_fasta.sh {input.pangraph_fasta} $(cat {input.gene_block}) > {output.gene_block_fasta}")
		shell("makeblastdb -in {output.gene_block_fasta} -dbtype 'nucl'")
		shell("blastn -max_hsps 10000 -query {input.gene_fasta} -db {output.gene_block_fasta} -outfmt '6 sstart send slen' > {output.gene_rel_locations}")
		shell("python scripts/compute_distances.py --block_csv {input.block_csv} \
											--gene_block_file {input.gene_block} \
											--snps {input.snps}\
											--output {output.dists}\
											--gene_offset {output.gene_rel_locations}")

rule positional_entropies:
	input:
		pangraph="output/pangraph/{gene}/{gene}_pangraph.json",
	params:
		blastdb_fasta="output/pangraph/{gene}/{gene}_extracted.fa",
		gene_query="input/focal_genes/{gene}.fa",
		gene_locations="output/pangraph/{gene}/{gene}_gene_locations.txt",
		assignments="output/gene_variants/sequence_assignments/{gene}.csv"
	output:
		real="output/pangraph/{gene}/positional_entropies.txt",
		consensus="output/pangraph/{gene}/positional_entropies_consensus.txt",
		consensus_relative="output/pangraph/{gene}/positional_entropies_consensus_relative.txt"
	run:
		shell("makeblastdb -in {params.blastdb_fasta} -dbtype 'nucl'")
		shell("blastn -max_hsps 10000 -query {params.gene_query} -db {params.blastdb_fasta} -outfmt '6 sseqid sstart send' > {params.gene_locations}")
		shell("python scripts/positional_entropy.py --json {input.pangraph} --normalise > {output.real}")
		shell("python scripts/positional_entropy.py --json {input.pangraph} --normalise --consensus --genelocations {params.gene_locations} > {output.consensus_relative}")
		shell("python scripts/positional_entropy.py --json {input.pangraph} --normalise --consensus > {output.consensus}")

rule plot_breakpoint_distances:
	input:
		dists="output/pangraph/{gene}/{gene}.output_dists.csv",
		deduplicated_gene="output/gene_variants/sequence/{gene}_seqs.fa.dedup.txt",
		variant_assignments="output/gene_variants/sequence_assignments/{gene}.csv"
	params:
		gene="{gene}",
		output_pdf_prefix="output/pangraph/{gene}/plots/{gene}_breakpoint_distances"
	output:
		"output/pangraph/{gene}/plots/{gene}_breakpoint_distances-all.pdf",
	shell:
		"Rscript scripts/plot_output_dists.R {input.dists} {input.deduplicated_gene} {input.variant_assignments}\
						--output_pdf_prefix {params.output_pdf_prefix} \
						--focal_gene {params.gene}"

rule plot_linear_blocks:
	input:
		block_csv="output/pangraph/{gene}/{gene}_pangraph.json.blocks.csv",
		focal_block_file="output/pangraph/{gene}/{gene}.gene_block.txt"
	output:
		"output/pangraph/{gene}/plots/linear_blocks.pdf"
	run:
		shell("Rscript scripts/plot_blocks_linear.R {input.block_csv} --focal_block {input.focal_block_file} --output_pdf {output}")

rule plot_bandage:
	input:
		"output/pangraph/{gene}/{gene}_pangraph.gfa.coloured.gfa"
	params:
		output_png="output/pangraph/{gene}/plots/{gene}_pangraph.bandage_plot.png"
	output:
		"output/pangraph/{gene}/plots/bandage.log_file"
		#"output/pangraph/{gene}/plots/{gene}_pangraph.bandage_plot.png"
	shell:
		"Bandage image {input} {params.output_png} --height 4000 --width 7000 --colour custom > {output}" if config["bandage"]==True else
		"echo 'Bandage not run' > {output}"

# rule combine_linear_and_bandage:
# 	input:
# 		linear="output/pangraph/{gene}/plots/linear_blocks.pdf",
# 		bandage_logfile="output/pangraph/{gene}/plots/bandage.log_file"
# 	params:
# 		output_pdf = "output/pangraph/{gene}/plots/{gene}_bandage_and_linear.pdf",
# 		bandage="output/pangraph/{gene}/plots/{gene}_pangraph.bandage_plot.png"
# 	output:
# 		"output/pangraph/{gene}/plots/{gene}_bandage_and_linear.logfile"
# 	shell:
# 		"Rscript scripts/combine_two_plots.R {input.linear} {params.bandage} --output_pdf {params.output_pdf} > {output}" if config["bandage"]==True else
# 		"echo 'Bandage not run' > {output}"


rule plot_positional_entropies_consensus_relative:
	input:
		"output/pangraph/{gene}/positional_entropies_consensus_relative.txt"
	output:
		"output/pangraph/{gene}/plots/{gene}_positional_entropies_consensus_relative.pdf"
	run:
		shell("Rscript scripts/plot_entropies.R {input} --relative T --output_pdf {output}")

rule plot_NJ_tree_central_gene:
	input:
		tree="output/gene_variants/sequence/{gene}_seqs.fa.dedup.aln.refined.tre",
		aln="output/gene_variants/sequence/{gene}_seqs.fa.dedup.aln",
		variant_assignments="output/gene_variants/sequence_assignments/{gene}.csv",
		dup_names="output/gene_variants/sequence/{gene}_seqs.fa.dedup.txt"
	output:
		"output/pangraph/{gene}/plots/NJ_tree_central_gene.pdf"
	shell:
		"Rscript scripts/plot_tree_variants.R --aln {input.aln} \
											--tree {input.tree} \
											--variants {input.variant_assignments}\
											--dup_names {input.dup_names} \
											--output_pdf {output}"

rule plot_output_dists_altair:
	input:
		"output/pangraph/{gene}/{gene}.output_dists.csv"
	output:
		"output/pangraph/{gene}/plots/{gene}_ecdf.html"
	shell: 
		"python scripts/plot_output_dists_altair.py --dist_csv {input} --output_html {output}"


rule plot_linear_blocks_altair:
	input:
		"output/pangraph/{gene}/{gene}_pangraph.json.blocks.csv"
	params:
		gene_name="{gene}", 
		gff="input/gffs/{gene}_annotations.gff",
		flanking_width=max(config["region_upstream"], config["region_downstream"])
	output:
		"output/pangraph/{gene}/plots/{gene}_linear_blocks.html",
	shell: 
		"python scripts/plot_linear_blocks_altair.py --flanking_width {params.flanking_width} --block_csv {input} --gene_name {params.gene_name} --output {output}" if config["include_gff"]==False else
		"python scripts/plot_linear_blocks_altair.py --flanking_width {params.flanking_width} --block_csv {input} --gene_name {params.gene_name} --output {output} --gff_file {params.gff}"
		
rule plot_linear_blocks_altair_deduplicated:
	input:
		"output/pangraph/{gene}/{gene}_pangraph.json.blocks.csv"
	params:
		gene_name="{gene}", 
		gff="input/gffs/{gene}_annotations.gff",
		flanking_width=max(config["region_upstream"], config["region_downstream"])
	output:
		"output/pangraph/{gene}/plots/{gene}_linear_blocks_deduplicated.html"
	shell: 
		"python scripts/plot_linear_blocks_altair.py --flanking_width {params.flanking_width} --unique --block_csv {input} --gene_name {params.gene_name} --output {output}" if config["include_gff"]==False else
		"python scripts/plot_linear_blocks_altair.py --flanking_width {params.flanking_width} --unique --block_csv {input} --gene_name {params.gene_name} --output {output} --gff_file {params.gff}"

# rule combine_breakpoint_and_NJ:
# 	input:
# 		breakpoint_pdf="output/pangraph/{gene}/plots/plot_breakpoint_distances-all.pdf",
# 		NJ_tree="output/pangraph/{gene}/plots/NJ_tree_central_gene.pdf"
# 	params:
# 		output_pdf="output/pangraph/{gene}/plots/NJ_tree_breakpoint_combined.pdf"
# 	output:
# 		"output/pangraph/{gene}/plots/breakpoint_and_NJ.logfile"
# 	shell:
# 		"Rscript scripts/combine_two_plots.R {input.breakpoint_pdf} {input.NJ_tree} \
# 			--output_pdf {params.output_pdf} --rel_width 1.8 --width 10 --height 4 > {output}" if config["version"]=="default" else
# 		"echo 'Not run due to on cluster (no pdftools)' > {output}"
