
#
import itertools as iter
import argparse

def get_options():
    parser = argparse.ArgumentParser(description='computes a variety of distances between sequences')
    parser.add_argument("--block_csv", help="Input file (pangraph json)", type=str)
    parser.add_argument("--gene_block_file", help='File containing block with gene of interest (anchor)', type=str)
    parser.add_argument("--gene_offset", help='File containing actual locations of anchor gene in consensus central block', type=str)
    parser.add_argument("--snps", help="file of SNP distances between isolates", type=str)
    parser.add_argument("--gene_assignments", help="file of variant assignment of genes", type=str)
    parser.add_argument("--output", help="output file", type=str)
    return parser.parse_args()

def distFirstBreakpoint(a, b, starting_block, upstream=True):
    """a, b are lists of blocks"""
    a_start = a.index(starting_block)
    b_start = b.index(starting_block)
    shared_blocks = []
    if upstream==False:
        i = 0
        while i<len(a)-a_start-1 and i<len(b)-b_start-1:
            i += 1
            if b_start<len(b):
                if a[i+a_start]==b[b_start+i]:
                    shared_blocks.append(b[b_start+i])
                else:
                    #if a[i+a_start] in shared_blocks: # append the block if it's a duplicated pre-existing one
                    #    shared_blocks.append(a[i+a_start]) # this means duplications 'don't count'
                    break
    elif upstream==True:
        i = 0
        while a_start-i > 0 and b_start-i > 0:
            i = i+1
            if a[a_start-i]==b[b_start-i]:
                shared_blocks.append(b[b_start-i])
            else:
                #if a[a_start+i] in shared_blocks:
                #    shared_blocks.append(a[a_start+i])
                break
    return(shared_blocks)


def jaccard(a, b):
    """a, b are lists of blocks"""
    # Needs writing
    return


def main():
    args = get_options()

    gene_assignments = {line.strip().split(",")[0]:line.strip().split(",")[1] for line in open(args.gene_assignments, "r").readlines()}

    with open(args.gene_offset, 'r') as f:
        gene_offset = f.read().strip('\n').split('\t')
        # formatted as balst output of focal gene against the consensus gene block: sstart send slen 
    gene_offset_upstream = int(gene_offset[0])
    gene_offset_downstream = int(gene_offset[2])-int(gene_offset[1])
    block_lengths = {}
    path_dict = {}
    with open(args.block_csv, 'r') as f:
        for i, line in enumerate(f.readlines()):
            if i>0:
                line = line.split(',')
                if line[0] in path_dict.keys():
                    path_dict[line[0]].append(line[1])
                else:
                    path_dict[line[0]] = [line[1]]
                if line[1] in block_lengths.keys():
                    block_lengths[line[1]][line[0]] = abs(int(line[4])-int(line[3]))
                else:
                    block_lengths[line[1]] = {line[0]: abs(int(line[4])-int(line[3]))}
    # Generate snp-dists?
    starting_block = open(args.gene_block_file, "r").readline().strip("\n")

    with open(args.output, "w") as output_file:
        output_file.write('seq1,seq2,dist.up,dist.down,snps\n')
        with open(args.snps, 'r') as f: # generated with snp-dists -p from gene seqs
            for line in f.readlines():
                line = line.strip().split()
                a, b, snps = line[0], line[1], int(line[2])
                if a in path_dict.keys() and b in path_dict.keys(): # check if in pangraph (may not be if e.g. seq was too short)
                    upstream_blocks = distFirstBreakpoint(path_dict[a], path_dict[b], starting_block)
                    upstream_dist = sum([block_lengths[x][a] for x in upstream_blocks])
                    downstream_blocks = distFirstBreakpoint(path_dict[a], path_dict[b], starting_block, upstream=False)
                    downstream_dist = sum([block_lengths[x][a] for x in downstream_blocks])
                    output_file.write('%s,%s,%d,%d,%d,%s,%s\n' % (a, b, upstream_dist+gene_offset_upstream, downstream_dist+gene_offset_downstream, snps, gene_assignments[a], gene_assignments[b]))

if __name__ == "__main__":
    main()
