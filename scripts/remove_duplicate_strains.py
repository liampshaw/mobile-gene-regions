# given a list of strains and metadata, removes possibly duplicated strains
# i.e. those with the same year, country, genus
# example usage:
# python scripts/remove_duplicate_strains.py --metadata data/metadata.csv 
# 											--strain_list beta-lactamases/accs/VIM-1_accs.txt 
# 											--metadata_fields TaxGenus,Year

import pandas as pd
import argparse

def get_options():
	parser = argparse.ArgumentParser(description='Remove possibly duplicated strains to address sampling bias')
	parser.add_argument('--strain_list', help='Initial list of strains', type=str, required=True)   
	parser.add_argument('--metadata', help='Metadata file to use to subset', type=str, required=True)   
	parser.add_argument('--metadata_fields', help='Metadata fields to remove duplicates from (comma-separated)', type=str, required=True)   
	parser.add_argument('--output', help='output file', required=False, type=str, default='')
	return parser.parse_args()


def main():
	args = get_options()
	metadata_df = pd.read_csv(args.metadata, index_col=0)
	strain_list = [line.strip() for line in open(args.strain_list, "r").readlines()]
	strain_list = [s for s in strain_list if s in metadata_df.index]
	metadata_df_strains = metadata_df.loc[strain_list]

	metadata_fields = args.metadata_fields.split(",")


	metadata_df_filtered = metadata_df_strains[metadata_fields].drop_duplicates()
	if args.output=="":
		for l in list(metadata_df_filtered.index):
			print(l)
	else:
		with open(args.output, 'w') as f:
			for l in list(metadata_df_filtered.index):
				f.write("%s\n" % l)
	
if __name__=="__main__":
	main()
