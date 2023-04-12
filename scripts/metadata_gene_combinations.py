import pandas as pd

metadata = pd.read_csv(snakemake.input.metadata, index_col=0)

gene = snakemake.params.gene

#GENES = ["CMY",\
#			"CTX-M",\
#			"GES",\
#			"IMP",\
#			"KPC",\
#			"NDM",\
#			"PER",\
#			"SHV",\
#			"TEM",\
#			"VEB",\
#			"VIM"]


#for gene in GENES:
with open(snakemake.output[0], 'w') as f:
	hits = metadata.loc[[(str(gene) in x) for x in list(metadata["Gene.hits"])]].index
	for h in hits:
		f.write("%s\n" % h)
