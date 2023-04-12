## Data files

### Provided with repository

`accs.txt` - list of accessions for the n=7,225 contigs in the dataset. 

`BLDB-subset.tsv` - subset of data from the Beta-Lactamase DataBase (available [here](http://bldb.eu/)), manually compiled and reformatted. Currently there is no download option)  
 
`metadata.csv` - a metadata file for the contigs included in the dataset, including information such as associated BioSample and (where available) country, year, and taxonomy of the isolate. For details on how this was programmatically compiled, see the pipeline in `_metadata-preparation`. 

`focal_genes.fasta` - a fasta file of the 'focal genes' that we analyse the genomic contexts around (individual gene fastas are in `focal_genes` directory).
 
### Downloaded by snakemake pipeline 

Beta-lactamase variants are numerous and the nomenclature is confusing. The pipeline can be run with two databases of named variants which can differ very slightly from each other:

`CARD_db.fa` - CARD database of nucleotide sequences of beta-lactamase gene variants [source](https://card.mcmaster.ca/download/0/broadstreet-v3.2.6.tar.bz2).  
`NCBI_db.fa` - NCBI AMRFinderPlus database of nucleotide sequences of beta-lactamase gene variants [source](https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/3.11/2022-12-19.1/AMR_CDS).

These are both downloaded when the `snakemake` command is run from the parent directory. In practice, it doesn't make much difference, so the default is to use CARD. 

### Not provided or downloaded

To download fasta files for all contigs, something like the following:

```
mkdir fastadir
cd fastadir
while read acc;
do
	ncbi-acc-download --format fasta $acc
done < ../accs.txt
```

This `fastadir` should then be set in the `configs` file used.
