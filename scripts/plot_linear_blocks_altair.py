import altair as alt
import csv
import json
import pandas as pd
import re
import argparse
import numpy as np
from itertools import combinations

def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn
import scipy.cluster.hierarchy as sch
alt.renderers.enable('html')

def get_options():
    parser = argparse.ArgumentParser(description='Produce html plot of linear blocks')
    parser.add_argument('--block_csv', help='CSV of blocks, produced from pangraph pipeline.', type=str)
    parser.add_argument("--gff_file", help='GFF of annotations (can be concatenated files)', type=str, default='')
    parser.add_argument("--gene_name", help='Name of focal gene (matched with re from GFF)', type=str)
    parser.add_argument('--output', help='Output html (default: output.html)', type=str, default='output.html')
    parser.add_argument('--unique', help='Whether to plot unique block configurations (for a simpler/smaller plot)', default=False, action='store_true')
    parser.add_argument('--flanking_width', help='Size of flanking regions', type=int, default=5000)
    parser.add_argument('--strain_list', help='Only include comparisons between sequences named in this file', type=str, default='', required=False)   
    return parser.parse_args()

def calculate_jaccard_index(set1, set2):
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union

def pairwise_jaccard_index(dictionary):
    keys = list(dictionary.keys())
    n = len(keys)
    jaccard_matrix = [[0.0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        set1 = set(dictionary[keys[i]])
        set2 = set(dictionary[keys[j]])
        jaccard_index = calculate_jaccard_index(set1, set2)
        jaccard_matrix[i][j] = jaccard_index
        jaccard_matrix[j][i] = jaccard_index
    return jaccard_matrix, keys

def convert_annotations(input_df, focal_gene_search_name):
    combined_df = []
    for genome in set(input_df['seqid']):
      df = input_df[(input_df['seqid']==genome)]
      focal_product_hits = [x for x in list(df['product']) if focal_gene_search_name in x]
      if len(focal_product_hits)==0:
        pass
      elif len(focal_product_hits)>0:
        focal_product = focal_product_hits[0]
        if (len(focal_product)!=0):
          focal_gene_strand = list(df[df['product']==focal_product]['strand'])[0]
        if (focal_gene_strand=="+"):
          focal_gene_start = list(df[df['product']==focal_product]['start'])[0]
          focal_gene_end = list(df[df['product']==focal_product]['end'])[0]
          df = df.assign(new_start =df['start']-focal_gene_start,
                        new_end = df['end']-focal_gene_start,
                        new_strand = df['strand'])
        if (focal_gene_strand=="-"):
          # We are going to reverse-complement everything
          # Choose an arbitrary maximum length D
          # to reverse things from
          # 0...10---20........500 where gene on -ve strand 10-20
          # becomes after reverse-complementing
          # 0.................480+++490...500
          # i.e. we use a maximum length to offset and flip strands
          D = max(df['end'])
          strand_convert = {'+':'-', '-':'+'}
          df = df.assign(new_start = D-df['end'],
                        new_end = D-df['start'],
                        new_strand = [strand_convert[x] for x in list(df['strand'])])
          focal_gene_start = focal_gene_start = list(df[df['product']==focal_product]['new_start'])[0]
          df['new_start'] = df['new_start']-focal_gene_start
          df['new_end'] = df['new_end']-focal_gene_start
      if (len(combined_df)==0):
        combined_df = df
      else:
        combined_df = pd.concat([combined_df, df])
    return combined_df


def my_theme():
  return {
    'config': {
      'view': {'continuousHeight': 300, 'continuousWidth': 800},  # from the default theme
    }
  }


alt.themes.register('my_theme', my_theme)
alt.themes.enable('my_theme')

# Create the text explanation
#explanation = ["Blocks are coloured by inferred sequence homology (grey=unique block, only seen in one sequence).", "Click on a block to highlight that block (hold shift and click to select multiple blocks).",
#                "Double-click to unselect blocks.",
#             "Y-axis order of sequences is determined by pairwise Jaccard index in terms of blocks."]
#explanation_gff = ["Blocks are coloured by inferred sequence homology (grey=unique block, only seen in one sequence).", "Annotations of CDSs can be toggled on/off using checkbox at bottom.",
#                "Mouseover of block or CDS annotation gives more information."]
               

# Arguments
# blocks.csv
# gff file(s)
# gene name
def main():
  args = get_options()

  block_df = pd.read_csv(args.block_csv)

  # remove the information on genome location (might want to keep, just useful for gff annotation consistency)
  block_df['genome'] = [re.sub(':.*', '', x) for x in list(block_df['genome'])]
  # To do: if we keep these, can use them to give actual positions with mouseover
  #print(len(block_df))

  # Check for strain list
  if args.strain_list!='':
    seqs_to_include = list(pd.read_csv(args.strain_list, header=None)[0])
    block_df = block_df.loc[(block_df['genome'].isin(seqs_to_include))]
    #print(len(block_df))

  # Sort out the block colours for the plot
  block_colours = {k:v for k, v in zip(block_df['block'], block_df['colour'])}
  sorted_block_colour_dict = {k: v for k, v in sorted(block_colours.items(), key=lambda item: item[0])}
  # Create a list of colours in the sorted order of the blocks
  sorted_block_colours = list(sorted_block_colour_dict.values())


  # If unique plot requested, subset to unique genome block paths
  if args.unique==True:
    # definitely this is very slow and could be made faster (encode blocks as ints and use array)
    genome_block_paths = {x:','.join(block_df['block'][(block_df['genome']==x)]) for x in set(block_df['genome'])}
    # invert genome
    inverted_genome_block_paths = {}
    for k, v in genome_block_paths.items():
      if v in inverted_genome_block_paths.keys():
        inverted_genome_block_paths[v] += [k]
      else:
        inverted_genome_block_paths[v] = [k]
    # take first genome entry for each unique path as a representative
    unique_genome_reps = [genomes[0] for genomes in sorted(inverted_genome_block_paths.values())]
    with open(args.output+'.unique_genomes.txt', 'w') as f:
      print('writing unique genome representatives')
      for g in unique_genome_reps:
        f.write(g+'\n')
    unique_genome_rep_counts = {g:len(sorted(inverted_genome_block_paths.values())[i]) for i, g in enumerate(unique_genome_reps)}
    block_df = block_df.loc[[i for i in range(len(block_df)) if block_df['genome'][i] in unique_genome_reps]]
    # New genome names
    new_genome_name_dict = {g:g+' (n='+str(unique_genome_rep_counts[g])+')' for g in set(block_df['genome'])}
    #old_genome_names = 
    #ordered_genomes = 
    block_df['old_genome'] = block_df['genome']
    block_df['genome'] = [new_genome_name_dict[g] for g in block_df['genome']]
    # New order (how does it feel?)
    #block_df['order'] = pd.Categorical(block_df['genome'], categories=ordered_genomes, ordered=True)

  # Jaccard distance between sequences in terms of blocks
  # First, get the sets of blocks (don't care about order)
  block_sets = {x:set(block_df[block_df['genome']==x]['block']) for x in set(block_df['genome'])}
  # Then, get the distances
  jaccard_matrix, keys = pairwise_jaccard_index(block_sets)
  linkage_matrix = sch.linkage(jaccard_matrix, method='average')
  dendrogram = sch.dendrogram(linkage_matrix, labels=keys)
  ordered_genomes = [keys[i] for i in dendrogram['leaves']]
  # and use these to order the genomes
  block_df['order'] = pd.Categorical(block_df['genome'], categories=ordered_genomes, ordered=True)

  # Selection of block - aim is to click one and highlight all others on chart
  block_selection = alt.selection_point(fields=['block'], empty=True) 
  block_plot = alt.Chart(block_df).mark_bar().encode(
    x=alt.X('start:Q', title='Position (bp)', scale=alt.Scale(domain=[0, max(block_df['end'])])),
    x2=alt.X2('end:Q'),
    y=alt.Y('genome:O', title='Genome', sort=ordered_genomes),
    color=alt.Color('block:N', scale=alt.Scale(range=sorted_block_colours), legend=None),
    opacity=alt.condition(block_selection, alt.value(1), alt.value(0.2)),
    tooltip=['block:N', 'start:Q', 'end:Q', 'strand:N']
).add_params(
    block_selection,
).interactive()

  # If GFF, use it
  if args.gff_file!='':
    # Read in gff 
    gff_df = pd.read_csv(args.gff_file, comment='#', sep='\t', header=None)
    gff_df.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    gff_df_cds = gff_df[gff_df['type']=='CDS']

    gff_df_cds = gff_df_cds.assign(product=[re.sub(';.*', '', re.sub('.*product=', '', x)) for x in gff_df_cds['attributes']],
                                  name=[re.sub('ID=', '', re.sub(';.*', '', re.sub('.*Name=', '', x))) for x in gff_df_cds['attributes']])


    # Convert the annotations
    protein_name = re.sub("-.*", "", args.gene_name.upper()) # needed for matching to protein 
    # CTX-M-65 becomes CTX as search string - this could go wrong if including different genes. 
    # In testing phase!
    annotation_hits = convert_annotations(gff_df_cds, protein_name)

    # Convert plot positions to start 2000bp downstream (should be able to change this...)
    annotation_hits['new_start'] = annotation_hits['new_start']+args.flanking_width
    annotation_hits['new_end'] = annotation_hits['new_end']+args.flanking_width

    # Reduce to only those within plot limits of linear block plot
    limits = [-100, max(block_df['end'])]
    annotation_hits_filtered = annotation_hits[(annotation_hits['new_start']>limits[0]) & 
                                              (annotation_hits['new_end']<limits[1]*1.1)]
    # add a 'genome' variable for consistency with block_df plot 
    if args.unique==False:
      annotation_hits_filtered = annotation_hits_filtered.assign(genome=annotation_hits_filtered['seqid'])
      # Get rid of any not in block plot
      annotation_hits_filtered = annotation_hits_filtered.loc[[x for x in annotation_hits_filtered.index if annotation_hits_filtered['genome'][x] in list(set(block_df['genome']))]]
    elif args.unique==True:
      # hacky fix to make the names of annotation seqids match up with the 'GENOME (n=1)' renaming
      #print(list(annotation_hits_filtered.index))
      #print(annotation_hits_filtered.loc[325]['seqid'])
      annotation_hits_filtered = annotation_hits_filtered.loc[[i for i in list(annotation_hits_filtered.index) if annotation_hits_filtered.loc[i]['seqid'] in unique_genome_reps]]
      annotation_hits_filtered = annotation_hits_filtered.assign(genome=[new_genome_name_dict[g] for g in annotation_hits_filtered['seqid']])

    # Add a measure of transposase/integrase
    annotation_hits_filtered['annotation_type'] = np.where(annotation_hits_filtered['product'].str.contains(protein_name), protein_name, # hack to make it lexicographically first
                                  np.where(annotation_hits_filtered['product'].str.contains('transposase|integrase|recombinase'), 'Mobile (transposase/integrase)',
                                           np.where(annotation_hits_filtered['product'].str.contains('hypothetical'), 'Hypothetical protein',
                                            'Other protein')))
    #sorted(annotation_hits_filtered['annotation_type'])

    annotation_colours = {k:v for k, v in zip([protein_name, 
                                          'Mobile (transposase/integrase/recombinase)', 
                                          'Other protein', 
                                          'Hypothetical Protein'], 
                                          ['red', 'black', 'darkgrey', 'gray'])}
    sorted_annotation_colours_dict= {k: v for k, v in sorted(annotation_colours.items(), key=lambda item: item[0])}
    # Create a list of colours in the sorted order of the blocks
    sorted_annotation_colours = list(sorted_annotation_colours_dict.values())


    # these colours are custom order
    annotation_colours = ['gray', 'red', 'black', 'darkgrey' ] 
  
    # Selection of gene - aim is to click one and highlight all others on chart
    #CDS_selection = alt.selection_point(fields=['product'], empty=True)
    product_highlight = alt.selection_point(fields=['product'], on='mouseover', nearest=True)
    #product_highlight_opacity = alt.selection_single(bind=product_highlight, fields=['opacity'], init={'opacity': True})

   # Create a checkbox selection to toggle the CDS annotations from gff on/off
    input_checkbox_annotations = alt.binding_checkbox(name="CDS annotations (NCBI) ")
    checkbox_selection_annotations = alt.param(bind=input_checkbox_annotations)
    opacity_checkbox_condition_annotations = alt.condition(
        checkbox_selection_annotations,
        alt.value(1),
        alt.value(0)
    )

    gff_plot = alt.Chart(annotation_hits_filtered).mark_rule(size=6, clip=True).encode(
        x=alt.X('new_start:Q', title=''),
        x2=alt.X2('new_end:Q'),
        y=alt.Y('genome:O', title='Genome', sort=ordered_genomes),
        #opacity=alt.condition(CDS_selection, alt.value(1), alt.value(0.01)), # doesn't currently work
        tooltip=['product:N', 'start:Q', 'end:Q', 'strand:N', 'name:N'],
        color=alt.Color('annotation_type:O', scale=alt.Scale(range=sorted_annotation_colours), legend=alt.Legend(title='Annotation', columns=1, symbolLimit=0))
        #color=alt.condition(
        #product_highlight & checkbox_selection_annotations,
        #alt.value('yellow'),
        #alt.Color('annotation_type:O', scale=alt.Scale(range=sorted_annotation_colours), legend=alt.Legend(title='Annotation', columns=1, symbolLimit=0)))
    )#.add_selection(product_highlight)

    gff_plot.configure_legend(
        titleFontSize=20,  # Adjust the title font size as needed
        labelFontSize=18  # Adjust the label font size as needed
    )



    # Add arrowheads
    heads_forward = alt.Chart(annotation_hits_filtered[(annotation_hits_filtered['new_strand']=='+')]).mark_point(
        shape='triangle', 
        size=50, 
        angle=90
    ).encode(
        x=alt.X('new_end:Q', title=''),
        y=alt.Y('genome:O', title='Genome', sort=ordered_genomes),
        color=alt.Color('annotation_type:O', scale=alt.Scale(range=annotation_colours), legend=None),
        fill=alt.Color('annotation_type:O', scale=alt.Scale(range=annotation_colours), legend=None)

        #opacity=alt.condition(CDS_selection, alt.value(1), alt.value(0))
      )
    heads_reverse = alt.Chart(annotation_hits_filtered[(annotation_hits_filtered['new_strand']=='-')]).mark_point(
        shape='triangle', 
        size=50, 
        angle=270).encode(
        color=alt.Color('annotation_type:O', scale=alt.Scale(range=annotation_colours), legend=None),
            fill=alt.Color('annotation_type:O', scale=alt.Scale(range=annotation_colours), legend=None),
        x=alt.X('new_start:Q', title=''),
        y=alt.Y('genome:O', title='Genome', sort=ordered_genomes),
      )

    # Combine to get an arrow plot
    gff_plot = (gff_plot+heads_forward+heads_reverse).add_params(
        #CDS_selection
    ).interactive()

    
    gff_checkbox = gff_plot.add_params(
      checkbox_selection_annotations
      ).encode(
        opacity=opacity_checkbox_condition_annotations
      )
    combined_chart = alt.layer(block_plot, gff_checkbox).resolve_scale(color='independent')
  # If no gff, just linear block plot
  elif args.gff_file=='':
    combined_chart = block_plot

  # Save the plot
  combined_chart.save(args.output, format='html')
  

if __name__=="__main__":
  main()
