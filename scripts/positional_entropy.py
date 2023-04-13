# Calculate positional entropy along a pangraph away from a central gene
import json
import argparse
from random import seed
from random import randint
import convert_pangraph_to_block_list as cp
import numpy as np
from math import e

# --json ../output/GES-24/GES-24-mmseqs2-polish.all_u5000_d5000_pangraph.json
# --geneblock LXDNJCUWIP


def get_options():
    parser = argparse.ArgumentParser(description='Calculate positional entropy along a pangraph away from a central gene')
    parser.add_argument('--json', help='Input json file', type=str)
    parser.add_argument('--name', help='Prefix (if concatenating multiple files together)', type=str, required=False, default='')
    #parser.add_argument('--geneblock', help='Block to centre the analysis on', type=str, default='', required=False)
    parser.add_argument('--step', help='step size to calculate entropy at', required=False, default=10)
    parser.add_argument('--region_size', help='how far to go out away from gene', required=False, default=5000)
    parser.add_argument('--genelocations', help='file of blast hits for gene in seqs: id start end', default='', required=False)
    parser.add_argument('--normalise', help='whether to normalise entropy data by log(N)', action='store_true', default=False, required=False)
    parser.add_argument('--subset', help='(optional) file of strains to subset too', type=str, required=False, default='')
    parser.add_argument('--consensus', help='use consensus lengths of blocks, not real lengths', required=False, default=False, action='store_true')
    return parser.parse_args()

def shannonEntropy(labels, base=None, normalise=False):
  value,counts = np.unique(labels, return_counts=True)
  norm_counts = counts / counts.sum()
  base = e if base is None else base
  if normalise==True:
    return -(((norm_counts * np.log(norm_counts)/np.log(base)).sum())/np.log(len(labels)))
  else:
    return -(norm_counts * np.log(norm_counts)/np.log(base)).sum()

def main():
    args = get_options()
    json_file = str(args.json)
    with open(json_file, 'r') as f:
      pangraph_json = json.load(f)

    genome_dict_original = cp.extractPaths(pangraph_json['paths'])
    # Don't care about exact lengths of blocks, will just use the consensus (average?) length
    block_dict = {x['id']:len(x['sequence']) for x in pangraph_json['blocks']}
    block_nums = {x['id']:i for i, x in enumerate(pangraph_json['blocks'])}
    block_dict_num = {block_nums[x['id']]:len(x['sequence']) for x in pangraph_json['blocks']}

 
    possible_paths = [[block_nums[x[0]] for x in genome_dict_original[y]] for y in genome_dict_original.keys()]

    # If we have a subset, use it
    if args.subset!='':
        subset_strains = [line.strip('\n') for line in open(args.subset, 'r').readlines()]
        genome_dict = {k:genome_dict_original[k] for k in subset_strains if k in genome_dict_original.keys()}
    else:
        genome_dict = genome_dict_original
    #print(genome_dict.keys())
    #print(genome_dict)

    
    #print(genomes_as_int[0])


    # # If we have the gene block, then we can use this to pin upstream/downstream
    # if args.geneblock!='':
    #     gene_block_num = block_nums[args.geneblock]
    #     #print('gene block is:', gene_block_num)
    #     # use the gene block position to get the starting point in each genome
    #     starting_point_upstream = [block[2] for block in genome_dict[g] for g in genome_dict.keys() if block[0]==args.geneblock]
    #     starting_point_downstream = [block[3] for block in genome_dict[g] for g in genome_dict.keys() if block[0]==args.geneblock]

    #     #print(starting_point_upstream)
    #     #print(starting_point_downstream)
    #     for i in range(0, 5000, 50):
    #         if starting_point_upstream[0]-i<0:
    #             pass
    #         else:
    #             upstream_block_vector = [g[starting_point_upstream[j]-i] for j, g in enumerate(genomes_as_int)]# this is the vector we want
    #             print(args.name, 'upstream', i, shannonEntropy(upstream_block_vector, args.normalise))
    #     for i in range(0, 5000, 50):
    #         if starting_point_downstream[0]+i>max(len(g) for g in genomes_as_int):
    #             pass
    #         else:
    #             downstream_block_vector = [g[starting_point_downstream[j]+i] for j, g in enumerate(genomes_as_int)]# this is the vector we want
    #             print(args.name, 'downstream', i, shannonEntropy(downstream_block_vector, args.normalise))
            
    if args.genelocations!='':
        gene_locations = {}
        for line in open(args.genelocations, 'r').readlines():
            line = line.strip().split('\t')
            gene_locations[line[0]] = [int(line[1]), int(line[2])]

        # Check if all of the genome_dict keys are in the gene_locations
        # They aren't if the total contig is shorter than the length or they had multiple 
        # for g in genome_dict.keys():
        #     if not g in gene_locations.keys():
        #         print("WARNING: genome ", g, " is not in gene locations")
        genome_dict = {g:genome_dict[g] for g in genome_dict.keys() if g in gene_locations.keys()}
        # Use actual positions in genomes
        genomes_as_int = [None]*len(genome_dict.keys())
        if args.consensus==False:
            genome_i = 0
            for g, path in genome_dict.items():
                #print(g)
                #genome_as_int = [None] * max([p[3] for p in path])

                genome_as_int = [item for sublist in [[block_nums[p[0]]]*(p[3]-p[2]) for p in path ] for item in sublist]

                genomes_as_int[genome_i] = genome_as_int
                genome_i += 1
        elif args.consensus==True:
            genome_i = 0
            for g, path in genome_dict.items():
                #print(g)
                #print([p[0] for p in path])
    #            genome_as_int = [None] * max([p[3] for p in path])
                #genome_as_int = [None] * [block_dict_num[p[0]] for p in path]

                genome_as_int = [item for sublist in [[block_nums[p[0]]]*block_dict[p[0]] for p in path ] for item in sublist]

                genomes_as_int[genome_i] = genome_as_int
                genome_i += 1


        #print('gene block is:', gene_block_num)
        # use the gene block position to get the starting point in each genome
        starting_point_upstream = [gene_locations[g][0] for  g in genome_dict.keys()]
        starting_point_downstream = [gene_locations[g][1] for  g in genome_dict.keys()]
        #print(starting_point_upstream)
        #print(starting_point_downstream)
        for i in range(0, int(args.region_size), int(args.step)):
            if starting_point_upstream[0]-i<0:
                pass
            else:
                upstream_block_vector = [g[starting_point_upstream[j]-i] for j, g in enumerate(genomes_as_int)]# this is the vector we want
                if args.name!="":
                    output_string = args.name+" "
                else:
                    output_string = ""
                output_string = output_string+'upstream '+str(i)+' '+str(shannonEntropy(upstream_block_vector, normalise=args.normalise))
                print(output_string)
        for i in range(0, int(args.region_size), int(args.step)):
            if max([x+i for x in starting_point_downstream])>=min(len(g) for g in genomes_as_int):
                pass
            else:
                downstream_block_vector = [g[starting_point_downstream[j]+i] for j, g in enumerate(genomes_as_int)]# this is the vector we want
                if args.name!="":
                    output_string = args.name+" "
                else:
                    output_string = ""
                output_string = output_string+'downstream '+str(i)+' '+str(shannonEntropy(downstream_block_vector, normalise=args.normalise))
                print(output_string)
            
    else: # otherwise, just print absolute positions
        genomes_as_int = [None]*len(genome_dict.keys())
        if args.consensus==False:
            genome_i = 0
            for g, path in genome_dict.items():
                #print(g)
                #genome_as_int min= [None] * max([p[3] for p in path])

                genome_as_int = [item for sublist in [[block_nums[p[0]]]*(p[3]-p[2]) for p in path ] for item in sublist]

                genomes_as_int[genome_i] = genome_as_int
                genome_i += 1
        elif args.consensus==True:
            genome_i = 0
            for g, path in genome_dict.items():
                #print(g)
                #print([p[0] for p in path])
    #            genome_as_int = [None] * max([p[3] for p in path])
                #genome_as_int = [None] * [block_dict_num[p[0]] for p in path]
                genome_as_int = [item for sublist in [[block_nums[p[0]]]*block_dict[p[0]] for p in path ] for item in sublist]

                genomes_as_int[genome_i] = genome_as_int
                genome_i += 1
        for i in range(0, int(args.region_size)*2, int(args.step)):
            if i<min([len(g) for g in genomes_as_int]):
                block_vector = [g[i] for g in genomes_as_int]# this is the vector we want
                if args.name=='':
                    print(i, shannonEntropy(block_vector, normalise=args.normalise))
                else:
                    print(args.name, i, shannonEntropy(block_vector, normalise=args.normalise))






if __name__=="__main__":
    main()
