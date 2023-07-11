import argparse
import pandas as pd

def get_options():
    parser = argparse.ArgumentParser(description='combines distances between sequences')
    parser.add_argument("--breakpoint_distances", help='breakpoint distances from compute_breakpoint_distance_blocks.py', type=str)
    parser.add_argument("--gene_locations", help='locations of focal gene in sequences', type=str, required=False, default='')
    parser.add_argument("--minimap2_dists", help="minimap2 distances")
    parser.add_argument("--output", help="output file", type=str)
    return parser.parse_args()


def main():
    args = get_options()
    breakpoints_df = pd.read_csv(args.breakpoint_distances)
    gene_locations_df = pd.read_csv(args.gene_locations, header=None, sep='\t')
    gene_locations_df.columns = ["seq1", "start", "end"]
    merged_df = pd.merge(breakpoints_df, gene_locations_df, on='seq1')
    gene_locations_df["seq2"] = gene_locations_df["seq1"]
    merged_df = pd.merge(merged_df, gene_locations_df[["seq2", "start", "end"]], on='seq2')
    merged_df.columns = ["seq1", "seq2", "dist.up.1", "dist.up.2", 
                        "dist.down.1", "dist.down.2",
                        "start.gene.1", "end.gene.1",
                        "start.gene.2", "end.gene.2"]

    merged_df["start.1"] = merged_df["start.gene.1"]-merged_df["dist.up.1"]-1
    merged_df["end.1"] = merged_df["end.gene.1"]+merged_df["dist.down.1"]
    merged_df["start.2"] = merged_df["start.gene.2"]-merged_df["dist.up.2"]-1
    merged_df["end.2"] = merged_df["end.gene.2"]+merged_df["dist.down.2"]
    print(merged_df)

    minimap2_df = pd.read_csv(args.minimap2_dists, header=None, sep='\t')
    minimap2_df["seq1"] = minimap2_df[0]
    minimap2_df["seq2"] = minimap2_df[5]
    minimap2_df["start.1"] = minimap2_df[2]
    minimap2_df["end.1"] = minimap2_df[3]
    minimap2_df["start.2"] = minimap2_df[7]
    minimap2_df["end.2"] = minimap2_df[8]
    # merge
    minimap2_df['seq_combined'] = minimap2_df.apply(lambda row: tuple(sorted([row['seq1'], row['seq2']])), axis=1)
    merged_df['seq_combined'] = merged_df.apply(lambda row: tuple(sorted([row['seq1'], row['seq2']])), axis=1)

    # Merge the DataFrames on the new combined column
    merged_df = pd.merge(merged_df[["seq_combined", "seq1", "seq2", "start.1", "end.1", "start.2", "end.2"]], 
                        minimap2_df[["seq_combined", "seq1", "seq2", "start.1", "end.1", "start.2", "end.2" ]], on='seq_combined')
    merged_df.to_csv(args.output)

if __name__=="__main__":
    main()