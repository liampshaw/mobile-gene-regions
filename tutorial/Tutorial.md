# Tutorial on pipeline

This document aims to give a detailed explanation of how to use the pipeline to explore the neighbourhood of a particular gene. 

N.B. This was written before the helper `analyze-flanking-regions.py` script was developed. It applies to using the snakemake pipeline by writing the config file yourself. 

### Simple example

Example data for n=52 chromosomes and plasmids containing *bla*IMP-4 is included in `input`. This can be run with either of the following commands:

```
python analyze-flanking-regions.py \
        --contigs input/contigs/IMP-4_contigs.fa \
        --gene_fasta input/focal_genes/IMP-4.fa \
        --focal_gene_name IMP-4 \
        --output_dir output
```

Using default parameters, this will extract the +/- 5kb flanking regions, orient the contigs so *bla*IMP-4 is on the positive strand, run pangraph, and produce plots. It should take about one minute to run. 

Output will be placed in `output/IMP-4`.


### Running the pipeline for multiple genes

The snakemake pipeline enables analysis of multiple focal genes with the same parameters from a single command. 

For each focal gene, your two input files should be put in the following places (from the main directory):

* `input/focal_genes/{focal_gene}.fa`
* `input/contigs/{focal_gene}_contigs.fa`

You should then create a config file (see e.g. `configs/example_config.yaml` in the main repo) and specify 

```
FOCAL_GENES = ["gene1", "gene2", "gene3"]
``` 

This means the pipeline will assume that the following files exist for each gene:

* `input/focal_genes/gene1.fa`
* `input/contigs/gene1_contigs.fa`

The snakemake pipeline is split up into different rules. You can run the full pipeline with:

```
snakemake --cores 1 --configfile your_config.yaml -r prepare_DB run_pangraph calculate_distances make_plots
```

The `beta-lactamases` folder contains a version of the pipeline specifically for twelve-betalactamase genes, run with a single command like this using the `betalactamases_config.yaml` file there. 

## Generating your own input files

The first thing you need to do is obtain two input files: a fasta file containing a focal gene of interest (`{gene}.fa`), and a set of contigs that contain the gene (`{gene}_contigs.fa`). 

You may already have these prepared for your own dataset. If so, you can run the pipeline with the helper script `analyze-flanking-regions.py`. Or, if you have multiple files, you can do so with snakemake and editing the config file as specified above. 

Alternatively, you may be interested in downloading some existing fasta files from NCBI. Here, we walk through the process of how to download contigs containing a particular antimicrobial resistance (AMR) gene.  


### Using NCBI to download AMR genes

Let us assume that you have a particular resistance gene you are interested in and want to look at the region around that gene in published assemblies. 

#### MicroBIGGE

The NCBI Pathogen Detection Microbial Browser for Identification of Genetic and Genomic Elements ([MicroBIGGE](https://www.ncbi.nlm.nih.gov/pathogens/microbigge/)) offers an interface to search for contigs containing a gene. 

For example, we can search for all contigs that contain *mcr-1.1* (a variant of the *mcr-1* mobilised colistin resistance gene) using the search `element_symbol:mcr-1.1`:

<img src="images/mcr1.1_search.png"  width="800">


**A note on searching for gene variants:** unfortunately the naming schemes for AMR genes can make this kind of Boolean search difficult. Taking *mcr-1* as an example, its variants are named as *mcr-1.1*, *mcr-1.2* etc. These are all almost identical with only a few SNPs separating them. You might want to include in a search using `element_symbol:mcr-1.*`. However, some hits will be to sequences that do not match a named variant, so will just be stored as *mcr-1*, which will not be matched by this search. But `element_symbol:mcr-1*` is not what you want either, because this would also return *mcr-10* - which is non-homologous to *mcr-1* (50.9% protein similarity). In conclusion, the right search to return all *mcr-1* variants but exclude other *mcr* genes is: `element_symbol:mcr-1.* OR element_symbol:mcr-1` (8987 hits as of 25 May 2023).


As of 25 May 2023 there was only an option to download whole contigs for up to 1,000 hits. For *mcr-1.1* there are 8403 hits, so we cannot download all contigs from this interface. 

Options for download are: 

* download the table of hits and download the contigs using the `Contig` accessions (e.g. with [ncbi-acc-download](https://github.com/kblin/ncbi-acc-download): `ncbi-acc-download -F fasta {accession}`. 
* download the gene hits with up to 2,000 bases flanking them in either direction (shorter contigs will return the whole contig)

The nice thing about MicroBIGGE is that you get access to the species information and some other metadata such as BioSample, which you can use to fetch other metadata. 

#### MicroBIGGE on Google Cloud Platform

An alternative way to search NCBI and download whole contigs is via the Google Cloud Platform. NCBI have a [guide](https://www.ncbi.nlm.nih.gov/pathogens/docs/microbigge_gcp/) on how to do this. 

We can use an SQL query to obtain a list of Google Cloud storage locations of gzipped contigs:

```
SELECT contig_acc
FROM `ncbi-pathogen-detect.pdbrowser.microbigge`
WHERE element_symbol LIKE 'mcr-1.%' OR element_symbol='mcr-1'
```

(returns 8987 hits as of 25 May 2023 i.e. matches MicroBIGGE)

This will return a list of the accessions of contigs, which we can download with `ncbi-acc-download -F fasta {accession}`. These contigs should then be combined into a single multi-fasta contig file for the pipeline (`gene_contigs.fa`).

Alternatively, you can select `contig_url` to get a list of Google Cloud storage locations for the contigs:

```
SELECT contig_url
FROM `ncbi-pathogen-detect.pdbrowser.microbigge`
WHERE element_symbol LIKE 'mcr-1.%' OR element_symbol='mcr-1'
```

We can download that list of locations as a csv (save as `locations.csv`):

<img src="images/gcp_save_csv.png"  width="800">

Then, having [installed](https://cloud.google.com/storage/docs/gsutil_install) `gsutil` which allows download of Cloud Storage locations from the command line, we can download the contigs:

```
while read f;
do
	gsutil cp $f .
done < locations.csv
``` 

### Optional: providing GFF files

The pipeline runs without annotation files by default, because pangraph uses only sequence similarity and no annotation information. However, it can be useful and interesting to see how the pangraph blocks correspond to the annotation information for the final plot of linear blocks. Only `CDS` features are used. 

If you have your own annotation files then they should be provided as a single gff. 

If you want to use NCBI annotation files, then you can get them with e.g. `ncbi-acc-download -F gff3 {accession}`. These can then be combined into a single gff (`cat *.gff > all.gff`; no need to strip out the headers etc., the pipeline just ignores all of that). 

If you want gff annotations on top of the linear blocks, set `include_gff: True` and put your combined gff in `input/gffs/{gene}_annotations.gff`.

WARNING: a large GFF (particularly where you have whole chromosomes) will take a long time to process for this final plot.

You can run the pipeline iteratively to get a better handle on this: first, run with `include_gff: False`. Then, look in `output/pangraph/{gene}/plots/{gene}_linear_blocks_deduplicated.html.unique_genomes.txt` to see the (reduced) list of accessions you would need gffs for to cover the variation in your dataset. Having downloaded these, you can then plot the reduced dataset with gff annotations with:

```
python scripts/plot_linear_blocks_altair.py --flanking_width {flanking_width} \
                                            --block_csv output/pangraph/{gene}/{gene}_pangraph.json.blocks.csv \
                                            --strain_list {deduplicated_genomes}
                                            --gene_name {gene} \ 
                                            --gff_file {your_combined_gff} \
                                            --output gff_plot.html
```


### Outputs

All output files are put into the `output_dir/{gene}` folder. 

Outputs to do with extracting the variation in the focal gene are put in `output_dir/{gene}/gene_diversity`. This is mainly relevant for SNV comparisons to the focal gene (TO DO: add integration with CARD or custom database.).

Pangraph outputs in `output_dir/{gene}/pangraph/`:
* `{gene}_pangraph.json` - pangraph of the flanking regions in all sequences
* `{gene}_pangraph.gfa` - gfa of pangraph. Can be inspected in Bandage
* `{gene}_pangraph.json.blocks.csv` - csv of block start/end positions in each sequence with colours assigned to show homology
* `{gene}.output_dists.csv` - breakpoint distances upstream/downstream between every pair of sequences
* `positional_entropies.txt` - block positional entropy at each position

Plots are outputted in `output_dir/{gene}/plots`:
* `NJ_tree_central_gene.pdf` - NJ tree of the focal gene
* `{gene}_pangraph.bandage_plot.png` - Bandage plot of pangraph
* `linear_blocks.pdf` - linear representation of pangraph 
* `{gene}_breakpoint_distances-...pdf` - empirical cumulative distribution function of pairwise breakpoint distances, either between `all`, `compare-to-focal-gene` (only comparisons involving focal gene), `compare-to-most-common` (only comparisons involving most common gene variant).
* `{gene}_linear_blocks.html` - interactive html of linear representation of pangraph (also `{gene}_linear_blocks_deduplicated.html` showing only unique patterns of blocks, with number of sequences having that pattern)
* `{gene}_ecdf.html` - interative html of all pairwise breakpoint distances


![](images/ecdf_screenshot.png)

![](images/linear_blocks_screenshot.png)

*To add*: more information about how to interact with linear blocks plot.










