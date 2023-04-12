# Metadata processing

After downloading raw metadata from NCBI for contigs, we process it to get a clean metadata file. 

Process involves (todo: make this as a `snakemake` pipeline):

* `01_download-genbank.sh`: downloads genbank data for list of NCBI nucleotide accessions
* `02_get-metadata.py`: downloads metadata for linked NCBI biosamples with `epost` (xml forat)
* `03_extract-xml.py`: extracts relevant xml into csv
* `04_combine-metadata.py`: combines csvs together to produce `metadata-raw.csv`
* `05a_process-metadata.R`: `metadata-raw.csv` > `metadata-processed-v1.csv`
* Manually edit spreadsheet: `metadata-processed-v1.csv` > `metadata-processed-v2.csv`
* `05b_process-metadata.R`: `metadata-processed-v2.csv` > `metadata-processed-v3.csv`

That final `v3` file becomes `data/metadata.csv` 

N.B. I've recently discovered [MicroBIGG-E](https://www.ncbi.nlm.nih.gov/pathogens/docs/microbigge_gcp/) which offers a way to download this data using SQL queries through Google BigQuery. I have not yet used this. 
