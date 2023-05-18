import altair as alt
import csv
import json
import pandas as pd

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

block_plot = alt.Chart(block_df).mark_bar().encode(
    x=alt.X('start:Q', title='Position'),
    x2=alt.X2('end:Q'),
    y=alt.Y('genome:O', title='Genome', sort=ordered_genomes),
    color=alt.Color('block:N', scale=alt.Scale(range=sorted_colours), legend=None),
    opacity=alt.condition(block_selection, alt.value(1), alt.value(0.3)),
    tooltip=['block:N', 'start:Q', 'end:Q', 'strand:N']
).add_params(
    block_selection
).interactive()#.properties(
 #   width=800,
#    height=200
#)

block_plot.save('blocks.html')



# Read in gff 
gff_df = pd.read_csv('NZ_PJHN01000019.gff', comment='#', sep='\t', header=None)
gff_df.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
gff_df_cds = gff_df[gff_df['type']=='CDS']

gff_df_cds = gff_df_cds.assign(product=[re.sub(';.*', '', re.sub('.*product=', '', x)) for x in gff_df_cds['attributes']],
                              name=[re.sub('ID=', '', re.sub(';.*', '', re.sub('.*Name=', '', x))) for x in gff_df_cds['attributes']])

# Convert the annotations
# Convert to positive strand for the focal gene
def convert_annotations(input_df, focal_gene):
  combined_df = []
  for genome in set(input_df['seqid']):
    df = input_df[input_df['seqid']==genome]
    focal_product = [x for x in list(df['product']) if 'MCR' in x][0]
    if (len(focal_product)!=0):
      focal_gene_strand = list(df[df['product']==focal_product]['strand'])[0]
    if (focal_gene_strand=="+"):
      focal_gene_start = list(df[df['product']==focal_product]['start'])[0]
      focal_gene_end = list(df[df['product']==focal_product]['end'])[0]
      df = df.assign(new_start =df['start']-focal_gene_start,
                    new_end = df['end']-focal_gene_start,
                    new_strand = df['strand'])
  return df

annotation_hits = convert_annotations(gff_df_cds, 'MCR')

# NEED TO ADD THE HANDLING OF NEGATIVE CASE BELOW    
  #   if (focal.gene.strand=="-"):
  #     # We are going to reverse-complement everything
  #     # Choose an arbitrary maximum length D
  #     # to reverse things from
  #     # 0...10---20........500 where gene on -ve strand 10-20
  #     # becomes after reverse-complementing
  #     # 0.................480+++490...500
  #     # i.e. we use a maximum length to offset and flip strands
  #     D = max(df$V5)
  #     df$new.start = D-df$V5
  #     df$new.end = D-df$V4
  #     df$new.strand = ifelse(df$V7=="-", "+", "-")
  #     focal.gene.start = df$new.start[grep(focal_gene, df$product)]
  #     df$new.start = df$new.start-focal.gene.start
  #     df$new.end = df$new.end-focal.gene.start
    
  #   df$forward = ifelse(df$new.strand=="+", T, F)
  #   if (is.null(combined_df)):
  #     combined_df = df
  #   else:
  #     combined_df = rbind(combined_df, df)
    
  # return(combined_df)

# 

# Reduce to only those with annotations provided
#annotation.hits = annotation.hits[which(annotation.hits$V1 %in% GENOMES),]

#genome.blocks = genome.blocks[which(genome.blocks$genome.ordered %in% GENOMES),]
#genome.blocks$genome.ordered =ordered(genome.blocks$genome.ordered, )

limits = [-max(block_df['end']), max(block_df['end'])]

annotation_hits_filtered = annotation_hits[(annotation_hits['new_start']>limits[0]) & 
                                          (annotation_hits['new_end']<limits[1]*5)]



gff_plot = alt.Chart(annotation_hits_filtered).mark_bar().encode(
    x=alt.X('start:Q', title='Position'),
    x2=alt.X2('end:Q'),
    y=alt.Y('seqid:O', title='Genome', sort=ordered_genomes),
    tooltip=['product:N', 'start:Q', 'end:Q', 'strand:N']
)

gff_plot.save('gff.html')
