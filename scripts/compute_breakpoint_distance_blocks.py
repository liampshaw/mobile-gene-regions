
#
import itertools as iter
import argparse

def get_options():
    parser = argparse.ArgumentParser(description='computes block distances between sequences in pangraph')
    parser.add_argument("--block_csv", help="Input file (produced from pangraph json)", type=str)
    parser.add_argument("--gene_block_file", help='File containing block with gene of interest (anchor)', type=str)
    parser.add_argument("--gene_offset", help='File containing actual locations of anchor gene in central block (needed for offset)', type=str, required=False, default='')
    parser.add_argument("--output", help="output file", type=str)
    return parser.parse_args()

def blocksFirstBreakpoint(a, b, starting_block, upstream=True):
    """a, b are lists of blocks. returns the blocks up until the first breakpoint"""
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



def main():
    args = get_options()

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
    # Block to start from (typically the one containing focal gene)
    starting_block = open(args.gene_block_file, "r").readline().strip("\n")

    # if 
    if args.gene_offset=='':
        with open(args.output, "w") as output_file:
            output_file.write('seq1,seq2,dist.up.1,dist.up.2,dist.down.1,dist.down.2\n')
            for a, b in iter.combinations(path_dict.keys(), 2):
                upstream_blocks = blocksFirstBreakpoint(path_dict[a], path_dict[b], starting_block)
                upstream_dist_a = sum([block_lengths[x][a] for x in upstream_blocks])
                upstream_dist_b = sum([block_lengths[x][b] for x in upstream_blocks])
                downstream_blocks = blocksFirstBreakpoint(path_dict[a], path_dict[b], starting_block, upstream=False)
                downstream_dist_a = sum([block_lengths[x][a] for x in downstream_blocks])
                downstream_dist_b = sum([block_lengths[x][b] for x in downstream_blocks])
                output_file.write('%s,%s,%d,%d,%d,%d\n' % (a, b, upstream_dist_a, upstream_dist_b, downstream_dist_a, downstream_dist_b))

    elif args.gene_offset!='':
        # if there is a gene offset given, use it
        gene_offsets = {}
        # assumes as blast output of focal gene against the extracted focal block:
        # sseqid sstart send  
        # so we add sstart to the upstream breakpoint distances 
        # and add (len_block-send) to the downstream breakpoint distancess
        with open(args.gene_offset, 'r') as f:
            for line in f.readlines():
                splitline = line.strip().split('\t')
                gene_offsets[splitline[0]] = [int(splitline[1]), block_lengths[starting_block][splitline[0]]-int(splitline[2])]
        with open(args.output, "w") as output_file:
            output_file.write('seq1,seq2,dist.up.1,dist.up.2,dist.down.1,dist.down.2\n')
            for seq1, seq2 in iter.combinations(path_dict.keys(), 2):
                # sort lexicographically
                a, b = sorted([seq1, seq2])
                a_offsets = gene_offsets[a]
                b_offsets = gene_offsets[b]
                upstream_blocks = blocksFirstBreakpoint(path_dict[a], path_dict[b], starting_block)
                upstream_dist_a = sum([block_lengths[x][a] for x in upstream_blocks])
                upstream_dist_b = sum([block_lengths[x][b] for x in upstream_blocks])
                downstream_blocks = blocksFirstBreakpoint(path_dict[a], path_dict[b], starting_block, upstream=False)
                downstream_dist_a = sum([block_lengths[x][a] for x in downstream_blocks])
                downstream_dist_b = sum([block_lengths[x][b] for x in downstream_blocks])
                output_file.write('%s,%s,%d,%d,%d,%d\n' % (a, b, upstream_dist_a+a_offsets[0], upstream_dist_b+b_offsets[0], downstream_dist_a+a_offsets[1], downstream_dist_b+b_offsets[1]))

    #gene_offset_upstream = int(gene_offset[0])
    #gene_offset_downstream = int(gene_offset[2])-int(gene_offset[1])
    


    
if __name__ == "__main__":
    main()
