import argparse
import pandas as pd


def get_options():
    parser = argparse.ArgumentParser(description='combine output distances for beta-lactamases')
    parser.add_argument('--block_dists',  help='file of output distances from blocks', required=True)
    parser.add_argument('--minimap2_dists',  help='file of output distances from minimap2', required=True)
    parser.add_argument('--output',  help='merged output csv', required=True)
    return parser.parse_args()

def main():
    args = get_options()
    block_dists = pd.read_csv(args.block_dists)
    #block_dists = block_dists.loc([y for y in [block_dists["seq1"][x]!=block_dists["seq2"][x] for x in range(0, len(block_dists))]])
    minimap2_dists = pd.read_csv(args.minimap2_dists, sep='\t', header=None)
    #print(block_dists)
    minimap2_dists = minimap2_dists[[0, 2, 3, 5]]
    minimap2_dists.columns = ["seq1", "start", "end", "seq2"]
    #print(minimap2_dists)

    minimap2_dists['seq_combined'] = minimap2_dists.apply(lambda row: tuple(sorted([row['seq1'], row['seq2']])), axis=1)
    block_dists['seq_combined'] = block_dists.apply(lambda row: tuple(sorted([row['seq1'], row['seq2']])), axis=1)

    # Merge the DataFrames on the new combined column
    merged_df = pd.merge(block_dists, minimap2_dists, on='seq_combined')

    # Remove the combined column from the merged DataFrame
    #merged_df.drop('seq_combined', axis=1, inplace=True)

    print(merged_df[["seq1_x", "seq2_x", "seq1_y", "seq2_y"]])
    print(merged_df[["dist.up", "dist.down", "start", "end"]])
    merged_df.to_csv(args.output)

if __name__=="__main__":
    main()