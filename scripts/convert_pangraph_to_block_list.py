# Convert pangraph json to list of blocks and positions
import json
import argparse
from random import seed
from random import randint



def get_options():
    parser = argparse.ArgumentParser(description='Convert pangraph json to path list of blocks with start/end positions')
    parser.add_argument('--json', help='Input json', type=str)
    parser.add_argument('--gfa', help='Input gfa', type=str) # TODO: remove gfa and generate with pangraph directly
    parser.add_argument('--all', '-a', help='colour all blocks (even unique ones)',
                    action='store_true')
    return parser.parse_args()


def strand(value): # converts True/False to +/- strand
    return ('' if value is None else
            '+' if value is True else
            '-' if value is False else
            str(value))


def colourBlocks(blocks, colour_all=False, colour_seed=1234):
    block_count_dict = {str(block):len(details['alignments']) for block,details in blocks.items()}
    # Create random colours
    block_colour_dict = {}
    if colour_all==False:
        block_num = len([x for x in block_count_dict.values() if x>1])
    if colour_all==True:
        block_num = len([x for x in block_count_dict.values()])
    block_colors = []
    seed(colour_seed) # for random colours
    for i in range(block_num):
        block_colors.append("#%06X" % randint(0, 0xFFFFFF))
    i = 0
    for block, count in block_count_dict.items():
        if count==1 and colour_all==False:
            block_colour_dict[block] = "#BDBABA" # colour grey
        else:
            block_colour_dict[block] = block_colors[i]
            i += 1
    return(block_colour_dict)

def extractPaths(extracted_json):
    genome_dict = {}
    paths = extracted_json['paths'] 
    path_id_to_name = {k: v["name"] for k,v in paths.items() }
    nodes = extracted_json['nodes']
    # to group nodes by path_id
    grouped = {}
    for node, node_details in nodes.items():
        #print(node, node_details)
        genome_name = path_id_to_name[str(node_details["path_id"])]
        block_info = [node_details["block_id"], node_details["strand"], node_details["position"][0], node_details["position"][1]]
        grouped.setdefault(genome_name, []).append(block_info)
    #print(grouped)
    return(grouped)
    #for path_entry in paths.items():
    #    path = path_entry[1] # to handle pangraph v1.2.0 json output
    #    print(path)
    #    genome = path['name']
    #    positions = path['position'] # N.B. already 1-indexed
    #    blocks = [x['id'] for x in path['blocks']]
    #    strands = [strand(x['strand']) for x in path['blocks']]
    #    genome_dict[genome] = []
    #    for i, b in enumerate(blocks):
    #        if i==len(blocks)-1:
    #            offset = 0 # From documentation: If N is the number of nodes this list has N+1 entries. The last entry is the position of the right edge of the last block in the path.
    #        else:
    #            offset = -1
    #        genome_dict[genome].append([b, strands[i], positions[i], positions[i+1]+offset])
    #return(genome_dict)

def rewriteGFA(gfa_file, colour_file, new_gfa_file):
    """Rewrites a GFA with colours."""
    block_colour_dict = {}
    with open(colour_file, 'r') as colours:
        for line in colours.readlines():
            block_colour_dict[line.split(',')[0]] = line.split(',')[1].strip('\n')
    with open(new_gfa_file, "w") as new_gfa:
        with open(gfa_file, 'r') as gfa:
            for line in gfa.readlines():
                entries = line.split()
                if entries[0] == "S":
                    entries.append("CL:z:"+block_colour_dict[entries[1]])
                new_gfa.write("%s\n" % "\t".join(entries))


def main():
    args = get_options()
    with open(str(args.json), 'r') as f:
        pangraph_json = json.load(f)

    output_blocks = str(args.json)+".blocks.csv"
    block_colour_file=str(args.json)+".colours.csv"
    # Get blocks and colour
    # see function above

    genome_dict = extractPaths(pangraph_json)
    block_colour_dict = colourBlocks(pangraph_json['blocks'])

    with open(output_blocks, 'w') as output_f:
        output_f.write("genome,block,strand,start,end,colour\n")
        for g, values in genome_dict.items():
            for v in values:
                block = str(v[0])
                output_f.write("%s\n" % (g+","+','.join([str(x) for x in v])+","+block_colour_dict[block]))
    with open(block_colour_file, 'w') as output_f_colours:
        output_f_colours.write("block,colour\n")
        for block, colour in block_colour_dict.items():
            output_f_colours.write("%s,%s\n" % (block, colour))

    rewriteGFA(str(args.gfa),
                block_colour_file,
                str(args.gfa)+".coloured.gfa")
    #rewriteGFA(gfa_file, colour_file, new_gfa_file):

    #print([ x for x in pangraph_json['paths']])
    #pangraph_json['paths'][0]['position']

if __name__=="__main__":
    main()
