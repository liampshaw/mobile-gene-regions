import altair as alt
import csv
import json
import pandas as pd
import re

from itertools import combinations
import scipy.cluster.hierarchy as sch


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


def my_theme():
  return {
    'config': {
      'view': {'continuousHeight': 300, 'continuousWidth': 800},  # from the default theme
    }
  }


alt.themes.register('my_theme', my_theme)
alt.themes.enable('my_theme')

block_df = pd.read_csv('blocks.csv')
# remove information on genome location (might want to keep, just useful for gff annotation consistency)
block_df['genome'] = [re.sub(':.*', '', x) for x in list(block_df['genome'])]

block_colours = {k:v for k, v in zip(block_df['block'], block_df['colour'])}
print(block_colours['LRQVYWPTXE'])
sorted_block_colours = {k: v for k, v in sorted(block_colours.items(), key=lambda item: item[0])}

# Create a list of colours in the sorted order of the blocks
sorted_colours = list(sorted_block_colours.values())


# Jaccard distance between sequences in terms of blocks
# First, get the sets of blocks (don't care about order)
block_sets = {x:set(block_df[block_df['genome']==x]['block']) for x in set(block_df['genome'])}
# Then, get the distances
jaccard_matrix, keys = pairwise_jaccard_index(block_sets)
linkage_matrix = sch.linkage(jaccard_matrix, method='average')
dendrogram = sch.dendrogram(linkage_matrix, labels=keys)
ordered_genomes = [keys[i] for i in dendrogram['leaves']]

block_df['order'] = pd.Categorical(block_df['genome'], categories=ordered_genomes, ordered=True)



# Selection of block - aim is to click one and highlight all others on chart
block_selection = alt.selection_point(fields=['block'], empty=True)





# Read in gff 
gff_df = pd.read_csv('mcr/all.gff', comment='#', sep='\t', header=None)
gff_df.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
gff_df_cds = gff_df[gff_df['type']=='CDS']

gff_df_cds = gff_df_cds.assign(product=[re.sub(';.*', '', re.sub('.*product=', '', x)) for x in gff_df_cds['attributes']],
                              name=[re.sub('ID=', '', re.sub(';.*', '', re.sub('.*Name=', '', x))) for x in gff_df_cds['attributes']])

# Convert the annotations
# Convert to positive strand for the focal gene
def convert_annotations(input_df, focal_gene):
  combined_df = []
  for genome in set(input_df['seqid']):
    df = input_df[(input_df['seqid']==genome)]
    focal_product = [x for x in list(df['product']) if 'MCR' in x][0]
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

annotation_hits = convert_annotations(gff_df_cds, 'MCR')

# Convert plot positions to start 2000bp downstream (should be able to change this...)
annotation_hits['new_start'] = annotation_hits['new_start']+2000
annotation_hits['new_end'] = annotation_hits['new_end']+2000

# Reduce to only those within plot limits of linear block plot
limits = [-1000, max(block_df['end'])]
annotation_hits_filtered = annotation_hits[(annotation_hits['new_start']>limits[0]) & 
                                          (annotation_hits['new_end']<limits[1]*1.1)]

annotation_hits_filtered = annotation_hits_filtered.assign(genome=annotation_hits_filtered['seqid'])

# Selection of gene - aim is to click one and highlight all others on chart
CDS_selection = alt.selection_point(fields=['product'], empty=True)

gff_plot = alt.Chart(annotation_hits_filtered).mark_rule(size=6, clip=True).encode(
    x=alt.X('new_start:Q', title='Position'),
    x2=alt.X2('new_end:Q'),
    y=alt.Y('genome:O', title='Genome', sort=ordered_genomes),
    #opacity=alt.condition(CDS_selection, alt.value(1), alt.value(0.01)), # doesn't currently work
    tooltip=['product:N', 'start:Q', 'end:Q', 'strand:N']
)

# Add arrowheads
heads_forward = alt.Chart(annotation_hits_filtered[(annotation_hits_filtered['new_strand']=='+')]).mark_point(
    shape='triangle', 
    size=50, 
    angle=90,
    fill='black',
    color='black').encode(
    x=alt.X('new_end:Q', title='Position'),
    y=alt.Y('genome:O', title='Genome', sort=ordered_genomes),
    #opacity=alt.condition(CDS_selection, alt.value(1), alt.value(0))
  )
heads_reverse = alt.Chart(annotation_hits_filtered[(annotation_hits_filtered['new_strand']=='-')]).mark_point(
    shape='triangle', 
    size=50, 
    angle=270,
       # opacity=alt.condition(CDS_selection, alt.value(1), alt.value(0)),
    fill='black',
    color='black').encode(
    x=alt.X('new_start:Q', title='Position'),
    y=alt.Y('genome:O', title='Genome', sort=ordered_genomes)
  )
# Combine to get an arrow plot
gff_plot = (gff_plot+heads_forward+heads_reverse).add_params(
    CDS_selection
).interactive()


gff_plot.save('gff.html')


#combined_chart = alt.layer(block_plot, gff_plot)

#combined_chart.save('gff_blocks.html')


# # Create a checkbox selection to toggle the CDS annotations from gff on/off
input_checkbox_annotations = alt.binding_checkbox(name="CDS annotations (NCBI) ")
checkbox_selection_annotations = alt.param(bind=input_checkbox_annotations)
opacity_checkbox_condition_annotations = alt.condition(
    checkbox_selection_annotations,
    alt.value(1),
    alt.value(0)
)

# For blocks - not currently working
# input_checkbox_blocks = alt.binding_checkbox(name="Make all blocks transparent ")
# checkbox_selection_blocks = alt.param(bind=input_checkbox_blocks)
# opacity_checkbox_condition_blocks = alt.condition(
#     checkbox_selection_blocks and block_selection,
#     alt.value(1),
#     alt.value(0.3)
# )


block_plot = alt.Chart(block_df).mark_bar().encode(
    x=alt.X('start:Q', title='Position (bp)', scale=alt.Scale(domain=[0, max(block_df['end'])])),
    x2=alt.X2('end:Q'),
    y=alt.Y('genome:O', title='Genome', sort=ordered_genomes),
    color=alt.Color('block:N', scale=alt.Scale(range=sorted_colours), legend=None),
    opacity=alt.condition(block_selection, alt.value(1), alt.value(0.2)),
    tooltip=['block:N', 'start:Q', 'end:Q', 'strand:N']
).add_params(
    block_selection
).interactive()#.properties(
 #   width=800,
#    height=200
#)

block_plot.save('blocks.html')

gff_checkbox = gff_plot.add_params(
    checkbox_selection_annotations
).encode(
    opacity=opacity_checkbox_condition_annotations
)
# Block checkbox - not currently working
#blocks_checkbox = block_plot.add_params(
#    checkbox_selection_blocks
#).encode(
#    opacity=opacity_checkbox_condition_blocks
#)

gff_checkbox.save('blocks_opacity.html')

combined_chart = alt.layer(blocks_plot, gff_checkbox)

# alt.concat(
#     combined_chart.encode(color='block:N'),
#     combined_chart.encode(color='product:N')
# ).resolve_scale(
#     color='independent'
# )
combined_chart.save('gff_blocks_checkbox.html')