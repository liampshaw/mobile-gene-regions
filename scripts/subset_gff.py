# For a formatted fasta file and a big pseudo-gff file,
# extract just the region of interest
import argparse
import pandas as pd
import re

def get_options():
	parser = argparse.ArgumentParser(description='Subset a big gff')
	parser.add_argument("--gff_file", help='GFF of annotations (can be concatenated files)', type=str, default='')
	parser.add_argument("--fasta_file", help='Fasta file (example header: >NZ_CP033439.1 10801bp 6694854-6705655 0 diffs)', type=str, default='')
	parser.add_argument("--padding", help='size of buffer zone (bases) to pad', type=str, default=1000)
	return parser.parse_args()

def main():
	args = get_options()
	gff_df = pd.read_csv(args.gff_file, comment='#', sep='\t', header=None)
	gff_df.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
	print(gff_df)
	genome_locations = {}
	with open(args.fasta_file, 'r') as f:
		for line in f.readlines():
			if line.startswith('>'):
				splitline = line.split(' ')
				name = re.sub(">", "", splitline[0])
				locations = [int(re.sub('[()]', '', re.sub('[a-z]', '', x))) for x in splitline[2].split('-')]
				genome_locations[name] = locations

	filtered_gff_df = pd.DataFrame(columns=gff_df.columns)
	# Iterate through each genome and its location
	for genome, locations in genome_locations.items():
		start_location, end_location = locations

		# Filter the gff_df DataFrame to get rows within the specified range
		filtered_rows = gff_df[(gff_df['seqid'] == genome) & 
	                           (gff_df['start']+args.padding >= start_location) &
	                           (gff_df['end']-+args.padding <= end_location)]

		# Append the filtered rows to the new DataFrame
		filtered_gff_df = pd.concat([filtered_gff_df, filtered_rows], ignore_index=True)

	print(len(gff_df))
	print(len(filtered_gff_df))
	print(filtered_gff_df)
	filtered_gff_df.to_csv("test.gff", sep="\t", header=None, index=False)




if __name__=="__main__":
	main()