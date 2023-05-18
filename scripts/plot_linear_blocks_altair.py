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



# Read in gff as well 