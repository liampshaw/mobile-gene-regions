# Produce alignment of a block
import argparse
import pandas as pd
import split_fasta as sf

def get_options():
    parser = argparse.ArgumentParser(description='extracts named block from sequences and aligns with mafft')
    parser.add_argument("--block_csv", help="Input file (produced from pangraph json)", type=str)
    parser.add_argument("--block_file", help='File containing name of block', type=str)
    parser.add_argument("--fasta", help="Input fasta that block_csv ultimately derived from")
    parser.add_argument("--output_fasta", help="output fasta")
    return parser.parse_args()


def main():
	args = get_options()
	block_name = [line.strip() for line in open(args.block_file, 'r').readlines()][0]

	blocks_df = pd.read_csv(args.block_csv)
	block_df = blocks_df[blocks_df['block']==block_name]
	#print(block_df)

	block_locations = {block_df['genome'].iloc[x] : [int(block_df['start'].iloc[x]), int(block_df['end'].iloc[x])] for x in range(0, len(block_df))}
	#print(block_locations)
	fasta = sf.read_fasta(args.fasta)
	with open(args.output_fasta, 'w') as f:
		for genome, block_location in block_locations.items():
			sequence = str(fasta[genome].seq)
			f.write('>%s %d:%d\n%s\n' % 
				(genome, block_location[0], block_location[1], sequence[(block_location[0]-1):(block_location[1]-1)]))


if __name__=="__main__":
	main()