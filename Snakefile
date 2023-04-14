# See configfile in e.g. configs/laptop_config.yaml for parameters
DB = config["DB"]
#FOCAL_GENE_DICT = config["focal_gene_dict"]
FOCAL_GENE_DICT = {"CMY-2": "CMY", 
					"CTX-M-15": "CTX-M",
					"CTX-M-65":"CTX-M",
					"GES-24":"GES",
					"IMP-4":"IMP",
					"KPC-2": "KPC",
					"NDM-1":"NDM",
					"PER-1":"PER",  
					"VIM-1":"VIM",
					"TEM-1": "TEM",
					"OXA-10": "OXA",
					"OXA-48":"OXA"}
rule all: # first rule - by default, runs pangraph (assumes have downloaded data etc.)
	input:
		expand("output/pangraph/{focal_gene}/{focal_gene}_pangraph.json", focal_gene=FOCAL_GENE_DICT.keys()),

rule gene_db:
	input:
		expand("DB/gene_alns/{gene}_{db}.tsv", gene=set(FOCAL_GENE_DICT.values()), db=DB),
		expand("DB/genes_plots/{db}-variants/{gene}.pdf", gene=set(FOCAL_GENE_DICT.values()), db=DB),


rule accessions:
	input:
		expand("output/analysis/{gene}_accessions.txt", gene=set(FOCAL_GENE_DICT.values()))

rule gene_db_assignments:
	input:
		expand("output/analysis/sequence_assignments/{focal_gene}_{db}.csv", focal_gene=FOCAL_GENE_DICT.keys(), db=DB)

rule seqs_extracted_snps:
	input:
		expand("output/analysis/sequence/{focal_gene}_seqs_extracted_from_contigs.snps.tsv", focal_gene=FOCAL_GENE_DICT.keys())

rule download_db:
	input:
		expand("data/{db}_db.fa", db=DB)

rule combined_contigs:
	input:
		expand("output/contigs/{gene}_combined_contigs.fa", gene=set(FOCAL_GENE_DICT.values())),
		expand("output/analysis/sequence/{focal_gene}_seqs_extracted_from_contigs.fa", focal_gene=FOCAL_GENE_DICT.keys()),
		expand("output/analysis/sequence/{focal_gene}_seqs_extracted_from_contigs.fa.dedup.aln", focal_gene=FOCAL_GENE_DICT.keys()),
		expand("output/analysis/sequence/{focal_gene}_seqs_extracted_from_contigs.fa.dedup.aln.refined.tre", focal_gene=FOCAL_GENE_DICT.keys()),

rule pangraph:
	input:
		expand("output/pangraph/{focal_gene}/{focal_gene}_extracted.fa", focal_gene=FOCAL_GENE_DICT.keys()),
		expand("output/pangraph/{focal_gene}/{focal_gene}_pangraph.gfa.coloured.gfa", focal_gene=FOCAL_GENE_DICT.keys()),
		expand("output/pangraph/{focal_gene}/{focal_gene}_{db}.output_dists.csv", focal_gene=FOCAL_GENE_DICT.keys(), db=DB),
		expand("output/pangraph/{focal_gene}/{focal_gene}_pangraph.json", focal_gene=FOCAL_GENE_DICT.keys()),
		expand("output/pangraph/{focal_gene}/{focal_gene}.gene_block.txt", focal_gene=FOCAL_GENE_DICT.keys())

rule plots:
	input:
		expand("output/pangraph/{focal_gene}/plots/{db}_plot_breakpoint_distances-all.pdf", focal_gene=FOCAL_GENE_DICT.keys(), db=DB),
		expand("output/pangraph/{focal_gene}/plots/bandage.log_file", focal_gene=FOCAL_GENE_DICT.keys()),
		expand("output/pangraph/{focal_gene}/plots/linear_blocks.pdf", focal_gene=FOCAL_GENE_DICT.keys()),
		expand("output/pangraph/{focal_gene}/plots/{focal_gene}_bandage_and_linear.logfile", focal_gene=FOCAL_GENE_DICT.keys()),
		#expand("output/pangraph/{focal_gene}/positional_entropies.txt", focal_gene=FOCAL_GENE_DICT.keys()),
		expand("output/pangraph/{focal_gene}/plots/positional_entropies.pdf", focal_gene=FOCAL_GENE_DICT.keys()),
		expand("output/pangraph/{focal_gene}/plots/positional_entropies_consensus.pdf", focal_gene=FOCAL_GENE_DICT.keys()),
		expand("output/pangraph/{focal_gene}/plots/positional_entropies_consensus_relative.pdf", focal_gene=FOCAL_GENE_DICT.keys()),
		expand("output/pangraph/{focal_gene}/plots/{db}_NJ_tree_central_gene.pdf", focal_gene=FOCAL_GENE_DICT.keys(), db=DB),
		expand("output/pangraph/{focal_gene}/plots/{db}_breakpoint_and_NJ.logfile", focal_gene=FOCAL_GENE_DICT.keys(), db=DB)



#############################
### SETTING UP DATABASES ####
#############################

rule download_NCBI_DB:
	output:
		"data/NCBI_db.fa"
	shell:
		"wget https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/3.11/2022-12-19.1/AMR_CDS -O {output}"

rule download_CARD_DB:
	output:
		"data/CARD_db.fa"
	run:
		shell("wget https://card.mcmaster.ca/download/0/broadstreet-v3.2.6.tar.bz2")
		shell("tar -xvf broadstreet-v3.2.6.tar.bz2 ./nucleotide_fasta_protein_homolog_model.fasta")
		shell("rm broadstreet-v3.2.6.tar.bz2")
		shell("mv nucleotide_fasta_protein_homolog_model.fasta {output}")

rule extract_genes_DB:
	input:
		fasta="data/{db}_db.fa"
	params:
		gene_name="{gene}",
		db="{db}"
	output:
		"DB/gene_fasta/{gene}_{db}.fa"
	shell:
		"python scripts/extract_gene_DB.py {input.fasta} {params.gene_name} {params.db} {output}"

rule align_genes_DB:
	input:
		"DB/gene_fasta/{gene}_{db}.fa"
	output:
		"DB/gene_alns/{gene}_{db}.aln"
	shell:
		"mafft --quiet --auto {input} > {output}"

rule snp_dists_DB:
	input:
		"DB/gene_alns/{gene}_{db}.aln"
	output:
		"DB/gene_alns/{gene}_{db}.tsv"
	shell:
		"snp-dists -a -b {input} > {output}"

rule run_fasttree_DB:
	input:
		"DB/gene_alns/{gene}_{db}.aln"
	output:
		"DB/gene_trees/{gene}_{db}.nwk"
	shell:
		"FastTree -quiet -nt -gtr {input} > {output}"

rule refine_fasttree_DB:
	input:
		tree= "DB/gene_trees/{gene}_{db}.nwk",
		aln= "DB/gene_alns/{gene}_{db}.aln"
	output:
		"DB/genes_refined-trees/{gene}_{db}.nwk"
	shell:
		"python scripts/refine_tree.py --tree_in {input.tree} --aln {input.aln} --tree_out {output}"

rule plot_tree_DB:
	input:
		"DB/genes_refined-trees/{gene}_{db}.nwk"
	output:
		"DB/genes_plots/{db}-variants/{gene}.pdf"
	shell:
		"Rscript scripts/plot_tree.R {input} {output}"

rule run_metadata:
	input:
		metadata="data/metadata.csv"
	params:
		gene="{gene}"
	output:
		"output/analysis/{gene}_accessions.txt"
	script:
		"scripts/metadata_gene_combinations.py"


rule copy_variant_fasta: # highly inelegant - copies files so is wasteful. But I don't understand symbolic links in snakemake...
    input:
        lambda wildcards: f"DB/gene_fasta/{FOCAL_GENE_DICT[wildcards.focal_gene]}"+"_{db}.fa"
    output:
        "DB/variant_fasta/{focal_gene}_{db}.fa"
    shell:
        "cp {input} {output}"


#############################
# GETTING CONTIG DATA READY #
#############################
rule combine_fastas_containing_gene:
	input:
		"output/analysis/{gene}_accessions.txt"
	output:
		"output/contigs/{gene}_combined_contigs.fa"
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
		variants="DB/variant_fasta/{focal_gene}_{db}.fa"
	output:
		"output/analysis/sequence_assignments/{focal_gene}_{db}.csv"
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
		"output/pangraph/{focal_gene}/{focal_gene}_{db}.output_dists.csv"
	shell:
		"python scripts/compute_distances.py --block_csv {input.block_csv} \
											--gene_block_file {input.gene_block} \
											--snps {input.snps}\
											--output {output}"

rule copied_assignments:
	input:
		expand("output/analysis/sequence_assignments/{focal_gene}.csv", focal_gene=FOCAL_GENE_DICT.keys())

rule copy_assignments:
	input:
		"output/analysis/sequence_assignments/{focal_gene}_CARD.csv"
	output:
		"output/analysis/sequence_assignments/{focal_gene}.csv"
	shell:
		"cp {input} {output}"

rule positional_entropies:
	input:
		pangraph="output/pangraph/{focal_gene}/{focal_gene}_pangraph.json",
	params:
		blastdb_fasta="output/pangraph/{focal_gene}/{focal_gene}_extracted.fa",
		gene_query="data/focal_genes/{focal_gene}.fa",
		gene_locations="output/pangraph/{focal_gene}/{focal_gene}_gene_locations.txt",
		assignments="output/analysis/sequence_assignments/{focal_gene}.csv"
	output:
		real="output/pangraph/{focal_gene}/positional_entropies.txt",
		consensus="output/pangraph/{focal_gene}/positional_entropies_consensus.txt",
		consensus_relative="output/pangraph/{focal_gene}/positional_entropies_consensus_relative.txt",
		consensus_relative_focal="output/pangraph/{focal_gene}/positional_entropies_consensus_relative_focal.txt",
		focal_subset="output/pangraph/{focal_gene}/{focal_gene}_focal_subset.txt"
	run:
		shell("makeblastdb -in {params.blastdb_fasta} -dbtype 'nucl'")
		shell("blastn -max_hsps 10000 -query {params.gene_query} -db {params.blastdb_fasta} -outfmt '6 sseqid sstart send' > {params.gene_locations}")
		shell("grep $(head -n 1 {params.gene_query} | tr -d '>') {params.assignments} | cut -d ',' -f 1 > {output.focal_subset}")
		shell("python scripts/positional_entropy.py --json {input.pangraph} --normalise > {output.real}")
		shell("python scripts/positional_entropy.py --json {input.pangraph} --normalise --consensus --genelocations {params.gene_locations} > {output.consensus_relative}")
		shell("python scripts/positional_entropy.py --json {input.pangraph} --subset {output.focal_subset} --normalise --consensus --genelocations {params.gene_locations} > {output.consensus_relative_focal}")
		shell("python scripts/positional_entropy.py --json {input.pangraph} --normalise --consensus > {output.consensus}")

#############################
##### PLOTTING OUTPUTS  #####
#############################

rule plot_positional_entropies:
	input:
		expand("output/pangraph/{focal_gene}/plots/positional_entropies.pdf", focal_gene=FOCAL_GENE_DICT.keys()),
		expand("output/pangraph/{focal_gene}/plots/positional_entropies_consensus.pdf", focal_gene=FOCAL_GENE_DICT.keys()),
		expand("output/pangraph/{focal_gene}/plots/positional_entropies_consensus_relative.pdf", focal_gene=FOCAL_GENE_DICT.keys()),
		expand("output/pangraph/{focal_gene}/plots/positional_entropies_consensus_relative_focal.pdf", focal_gene=FOCAL_GENE_DICT.keys())

rule plot_positional_entropies_real:
	input:
		"output/pangraph/{focal_gene}/positional_entropies.txt"
	output:
		"output/pangraph/{focal_gene}/plots/positional_entropies.pdf"
	run:
		shell("Rscript scripts/plot_entropies.R {input} --output_pdf {output}")

rule plot_positional_entropies_consensus:
	input:
		"output/pangraph/{focal_gene}/positional_entropies_consensus.txt"
	output:
		"output/pangraph/{focal_gene}/plots/positional_entropies_consensus.pdf"
	run:
		shell("Rscript scripts/plot_entropies.R {input} --output_pdf {output} ")

rule plot_positional_entropies_consensus_relative:
	input:
		"output/pangraph/{focal_gene}/positional_entropies_consensus_relative.txt"
	output:
		"output/pangraph/{focal_gene}/plots/positional_entropies_consensus_relative.pdf"
	run:
		shell("Rscript scripts/plot_entropies.R {input} --relative T --output_pdf {output}")

rule plot_positional_entropies_consensus_relative_focal:
	input:
		"output/pangraph/{focal_gene}/positional_entropies_consensus_relative_focal.txt"
	output:
		"output/pangraph/{focal_gene}/plots/positional_entropies_consensus_relative_focal.pdf"
	run:
		shell("Rscript scripts/plot_entropies.R {input} --relative T --output_pdf {output}")



rule breakpoint_distances_plot:
	input:
		expand("output/pangraph/{focal_gene}/plots/{db}_plot_breakpoint_distances-all.pdf", 
			focal_gene=FOCAL_GENE_DICT.keys(), db=DB)

rule plot_breakpoint_distances:
	input:
		dists="output/pangraph/{focal_gene}/{focal_gene}_{db}.output_dists.csv",
		deduplicated_focal_gene="output/analysis/sequence/{focal_gene}_seqs_extracted_from_contigs.fa.dedup.txt",
		variant_assignments="output/analysis/sequence_assignments/{focal_gene}_{db}.csv"
	params:
		focal_gene="{focal_gene}",
		output_pdf_prefix="output/pangraph/{focal_gene}/plots/{db}_plot_breakpoint_distances"
	output:
		"output/pangraph/{focal_gene}/plots/{db}_plot_breakpoint_distances-all.pdf",
	shell:
		"Rscript scripts/plot_output_dists.R {input.dists} {input.deduplicated_focal_gene} {input.variant_assignments}\
						--output_pdf_prefix {params.output_pdf_prefix} \
						--focal_gene {params.focal_gene}"

rule plot_linear_blocks:
	input:
		block_csv="output/pangraph/{focal_gene}/{focal_gene}_pangraph.json.blocks.csv",
		focal_block_file="output/pangraph/{focal_gene}/{focal_gene}.gene_block.txt"
	output:
		"output/pangraph/{focal_gene}/plots/linear_blocks.pdf"
	run:
		shell("Rscript scripts/plot_blocks_linear.R {input.block_csv} --focal_block {input.focal_block_file} --output_pdf {output}")

rule plot_bandage:
	input:
		"output/pangraph/{focal_gene}/{focal_gene}_pangraph.gfa.coloured.gfa"
	params:
		output_png="output/pangraph/{focal_gene}/plots/{focal_gene}_pangraph.bandage_plot.png"
	output:
		"output/pangraph/{focal_gene}/plots/bandage.log_file"
		#"output/pangraph/{focal_gene}/plots/{focal_gene}_pangraph.bandage_plot.png"
	shell:
		"Bandage image {input} {params.output_png} --height 4000 --width 7000 --colour custom > {output}" if config["bandage"]==True else
		"echo 'Bandage not run' > {output}"

rule combine_linear_and_bandage:
	input:
		linear="output/pangraph/{focal_gene}/plots/linear_blocks.pdf",
		bandage_logfile="output/pangraph/{focal_gene}/plots/bandage.log_file"
	params:
		output_pdf = "output/pangraph/{focal_gene}/plots/{focal_gene}_bandage_and_linear.pdf",
		bandage="output/pangraph/{focal_gene}/plots/{focal_gene}_pangraph.bandage_plot.png"
	output:
		"output/pangraph/{focal_gene}/plots/{focal_gene}_bandage_and_linear.logfile"
	shell:
		"Rscript scripts/combine_two_plots.R {input.linear} {params.bandage} --output_pdf {params.output_pdf} > {output}" if config["bandage"]==True else
		"echo 'Bandage not run' > {output}"

rule plot_NJ_tree_central_gene:
	input:
		tree="output/analysis/sequence/{focal_gene}_seqs_extracted_from_contigs.fa.dedup.aln.refined.tre",
		aln="output/analysis/sequence/{focal_gene}_seqs_extracted_from_contigs.fa.dedup.aln",
		variant_assignments="output/analysis/sequence_assignments/{focal_gene}_{db}.csv",
		dup_names="output/analysis/sequence/{focal_gene}_seqs_extracted_from_contigs.fa.dedup.txt"
	output:
		"output/pangraph/{focal_gene}/plots/{db}_NJ_tree_central_gene.pdf"
	shell:
		"Rscript scripts/plot_tree_variants.R --aln {input.aln} \
											--tree {input.tree} \
											--variants {input.variant_assignments}\
											--dup_names {input.dup_names} \
											--output_pdf {output}"


rule combine_breakpoint_and_NJ:
	input:
		breakpoint_pdf="output/pangraph/{focal_gene}/plots/{db}_plot_breakpoint_distances-all.pdf",
		NJ_tree="output/pangraph/{focal_gene}/plots/{db}_NJ_tree_central_gene.pdf"
	params:
		output_pdf="output/pangraph/{focal_gene}/plots/{db}_NJ_tree_breakpoint_combined.pdf"
	output:
		"output/pangraph/{focal_gene}/plots/{db}_breakpoint_and_NJ.logfile"
	shell:
		"Rscript scripts/combine_two_plots.R {input.breakpoint_pdf} {input.NJ_tree} \
			--output_pdf {params.output_pdf} --rel_width 1.8 --width 10 --height 4 > {output}" if config["version"]=="default" else
		"echo 'Not run due to on cluster (no pdftools)' > {output}"



