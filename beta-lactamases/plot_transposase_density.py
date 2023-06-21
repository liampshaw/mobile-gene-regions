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
    parser.add_argument('--dist_csv', help='CSV of output distances, produced from pangraph pipeline.', type=str)
    parser.add_argument("--gff_file", help='GFF of annotations (can be concatenated files)', type=str, default='')
    parser.add_argument("--gene_name", help='Name of focal gene (matched with re from GFF)', type=str)
    parser.add_argument('--output_html', help='Output html (default: output.html)', type=str, default='output.html')
    parser.add_argument('--unique', help='Whether to plot unique block configurations (for a simpler/smaller plot)', default=False, action='store_true')
    parser.add_argument('--flanking_width', help='Size of flanking regions', type=int, default=5000)
    parser.add_argument('--strain_list', help='Only include comparisons between sequences named in this file', type=str, default='', required=False)   
    return parser.parse_args()

def convert_csv_to_json(csv_input, json_output):
    with open(csv_input, newline='') as csvfile:
        # Read the CSV file into a dictionary
        reader = csv.DictReader(csvfile)
        # Convert the dictionary to a list of rows
        rows = list(reader)
    # Convert the list of rows to a JSON string
    json_str = json.dumps(rows)
    # Print the JSON string to file
    with open(json_output, 'w') as f:
            f.write(json_str)
    return
  
def map_to_categorical_snv(values):
    '''Maps SNVs to categories - have >7 as (arbitrary) highest category'''
    return_values = [""] * len(values)
    for i,v in enumerate(values):
        if v<=7:
            return_values[i] = str(v)
        else:
            return_values[i] = ">7"
    return return_values


def ecdf(values):
    '''returns the empirical cumulative distribution function for a list of values
    only keeps the steps in distribution to reduce size'''
    v_sorted = [0]+sorted(values)
    return pd.DataFrame({'ecdf':[x/len(values) for x in range(0, len(values)+1)],
                    'dist':v_sorted}).groupby('dist')['ecdf'].max().reset_index()



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
      # add relative start/end for upstream downstream
      #for _, row in df.iterrows():
      #  if row['new_start']<focal_gene_start:
      #df['relative_start'] = focal_gene_start-df['new_start']
      #df['relative_']

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
explanation = ["Blocks are coloured by inferred sequence homology (grey=unique block, only seen in one sequence).", "Click on a block to highlight that block (hold shift and click to select multiple blocks).",
                "Double-click to unselect blocks.",
             "Y-axis order of sequences is determined by pairwise Jaccard index in terms of blocks."]
explanation_gff = ["Blocks are coloured by inferred sequence homology (grey=unique block, only seen in one sequence).", "Annotations of CDSs can be toggled on/off using checkbox at bottom.",
                "Mouseover of block or CDS annotation gives more information."]
               

# Arguments
# blocks.csv
# gff file(s)
# gene name
def main():
  args = get_options()

  if args.gene_name!='':
    gene = args.gene_name.upper()
    df = pd.read_csv(args.dist_csv)
    df = df.loc[(df["variant1"]==gene) | (df["variant2"]==gene)]
    tmp_csv = args.dist_csv+'.'+gene+'.csv'
    df.to_csv(tmp_csv)
    # Convert csv to json
    tmp_json = tmp_csv+'_tmp.json'
    convert_csv_to_json(csv_input=tmp_csv, json_output=tmp_json)
    df = pd.read_json(tmp_json)

  else:
      # Convert csv to json
      tmp_json = args.dist_csv+'_tmp.json'
      convert_csv_to_json(csv_input=args.dist_csv, json_output=tmp_json)
      df = pd.read_json(tmp_json)

  # Check for strain list
  if args.strain_list!='':
      seqs_to_include = list(pd.read_csv(args.strain_list, header=None)[0])
      df = df.loc[(df['seq1'].isin(seqs_to_include)) & (df['seq2'].isin(seqs_to_include))]

  # Convert the data types of the columns to float 
  df['distup'] = df['dist.up'].astype(float)
  df['distdown'] = df['dist.down'].astype(float)
  # And add categorical SNVs
  df['snvscategorical'] = map_to_categorical_snv(list(df['snvs']))
  # Maximum distance
  # to clip the edges where ecdf artifically goes to 1 (so plot drops to 0) because we don't look beyond
  max_dist = max(df['distup'])-10 # 

  # For dropdown selection for SNV category of comparison
  options = sorted(set(df['snvscategorical']))
  labels = [option + ' ' for option in options]
  # Sort labels in the same order
  labels_sorted = [label for _, label in sorted(zip(options, labels))]

  # block_df = pd.read_csv(args.block_csv)

  # remove the information on genome location (might want to keep, just useful for gff annotation consistency)
  # block_df['genome'] = [re.sub(':.*', '', x) for x in list(block_df['genome'])]
  # To do: if we keep these, can use them to give actual positions with mouseover
   # Create ecdf dataframes
  plot_df_up = df.groupby('snvscategorical').apply(lambda x: ecdf(x['distup'])).reset_index(level=1, drop=True)
  plot_df_up = plot_df_up.assign(snvscategorical=plot_df_up.index)
  plot_df_down = df.groupby('snvscategorical').apply(lambda x: ecdf(x['distdown'])).reset_index(level=1, drop=True)
  plot_df_down = plot_df_down.assign(snvscategorical=plot_df_down.index)

  # If GFF, use it
  if args.gff_file!='':
    # Read in gff 
    gff_df = pd.read_csv(args.gff_file, comment='#', sep='\t', header=None)
    if args.strain_list!='':
      gff_df = gff_df[gff_df['seqid'].isin(strain_list)]
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

    # Reduce to only those within limits of output distances
    limits = [-100, 12000]
    annotation_hits_filtered = annotation_hits[(annotation_hits['new_start']>limits[0]) & 
                                              (annotation_hits['new_end']<limits[1]*1.1)]
    # add a 'genome' variable for consistency with block_df plot 
    if args.unique==False:
      annotation_hits_filtered = annotation_hits_filtered.assign(genome=annotation_hits_filtered['seqid'])
      # Get rid of any not in block plot
      #annotation_hits_filtered = annotation_hits_filtered.loc[[x for x in annotation_hits_filtered.index if annotation_hits_filtered['genome'][x] in list(set(block_df['genome']))]]
    elif args.unique==True:
      # hacky fix to make the names of annotation seqids match up with the 'GENOME (n=1)' renaming
      #print(list(annotation_hits_filtered.index))
      #print(annotation_hits_filtered.loc[325]['seqid'])
      annotation_hits_filtered = annotation_hits_filtered.loc[[i for i in list(annotation_hits_filtered.index) if annotation_hits_filtered.loc[i]['seqid'] in unique_genome_reps]]
      annotation_hits_filtered = annotation_hits_filtered.assign(genome=[new_genome_name_dict[g] for g in annotation_hits_filtered['seqid']])

    # Add a measure of transposase/integrase
    annotation_hits_filtered['annotation_type'] = np.where(annotation_hits_filtered['product'].str.contains(protein_name), protein_name, # hack to make it lexicographically first
                                  np.where(annotation_hits_filtered['product'].str.contains('transposase'), 'Transposase',
                                           np.where(annotation_hits_filtered['product'].str.contains('hypothetical'), 'Hypothetical protein',
                                            'Other protein')))

    # Determine the range of 'x' values based on the 'start' and 'end' columns
    min_x = annotation_hits_filtered['new_start'].min()
    max_x = annotation_hits_filtered['new_end'].max()
    all_x_values = range(int(min_x), int(max_x) + 1)

    # Initialize density dictionary with zeros for each 'annotation_type' and 'x' value
    density_dict = {k:{'x': [], 'density': []} for k in annotation_hits_filtered['annotation_type'].unique()}
    for k in density_dict.keys():
      for x in all_x_values:
        density_dict[k]['x'].append(x)
        density_dict[k]['density'].append(0)

    # Iterate over the rows of the original DataFrame and update density values
    for _, row in annotation_hits_filtered.iterrows():
        if row['annotation_type']=='Transposase':
          for x in range(int(row['new_start']), int(row['new_end']) + 1):
            density_dict['Transposase']['density'][x] += 1

    # Create the density_df DataFrame
    density_df = pd.DataFrame(density_dict['Transposase'])
    #density_df['density'] = density_df['density']-1 # 1 is default value so take off 1? Not sure about this...
    
    # at this point would like to convert to relative position
    # We have the gene positions
    median_gene_start = np.median(annotation_hits_filtered[annotation_hits_filtered['annotation_type']==protein_name]['new_start'])
    median_gene_end = np.median(annotation_hits_filtered[annotation_hits_filtered['annotation_type']==protein_name]['new_end'])
    # Define conditions and corresponding values
    conditions = [
        density_df['x'] <= median_gene_start,
        (density_df['x'] > median_gene_start) & (density_df['x'] <= median_gene_end),
        density_df['x'] > median_gene_start
    ]
    values = ['upstream', 'intermediate', 'downstream']
    # Assign location variable based on conditions
    density_df['location'] = np.select(conditions, values, default='unknown')

    #density_df = density_df.assign(location=['upstream' for x in list(density_df['x']) if x<median_gene_start else 'downstream']
    density_df = density_df[density_df['x']>0]
    density_df_up = density_df[density_df['x']<median_gene_start]
    density_df_down = density_df[density_df['x']>median_gene_end]
    density_df_up['dist'] = [5000-(x+1) for x in range(len(density_df_up))]
    density_df_down['dist'] = [x+1 for x in range(len(density_df_down))]

    # maximum density proportion
    max_density = max([density_df_up['density'].max(), density_df_down['density'].max()])
    # do 1 minus so that we can use the reversed axis of breakpoint distnace plot
    density_df_up['norm_density'] = [1-x/max_density for x in density_df_up['density']] 
    density_df_down['norm_density'] = [1-x/max_density for x in density_df_down['density']]


    # Set 'dist' as the index for just the 0 SNV comparisons
    plot_df_up_0 = plot_df_up[plot_df_up['snvscategorical']=="0"]
    plot_df_up_0 = plot_df_up_0.set_index('dist')
    # Create a new index with all integer values
    new_index = pd.Index(range(int(plot_df_up_0.index.min()), int(plot_df_up_0.index.max()) + 1), name='dist')
    # Reindex and fill in values
    plot_df_up_0 = plot_df_up_0.reindex(new_index, method='ffill').reset_index()

    # Ditto for downstream
    # Set 'dist' as the index for just the 0 SNV comparisons
    plot_df_down_0 = plot_df_down[plot_df_down['snvscategorical']=="0"]
    plot_df_down_0 = plot_df_down_0.set_index('dist')
    # Create a new index with all integer values
    new_index = pd.Index(range(int(plot_df_down_0.index.min()), int(plot_df_down_0.index.max()) + 1), name='dist')
    # Reindex and fill in values
    plot_df_down_0 = plot_df_down_0.reindex(new_index, method='ffill').reset_index()

    # Upstream
    ecdf_plot_up_unselected = alt.Chart(plot_df_up_0).mark_line(
        clip=True,
        strokeWidth=2,
        opacity=1,
    ).encode(
        x=alt.X("dist:Q",
                scale=alt.Scale(reverse=True, domain=[0, max_dist])),
        y=alt.Y('ecdf:Q', scale=alt.Scale(reverse=True, domain=[0,1])),
        color=alt.Color('snvscategorical:N', scale=None, 
        sort=alt.EncodingSortField('snvscategorical')),
        tooltip=[alt.Tooltip("dist:Q", title="")]
    )
    density_plot_up = alt.Chart(density_df_up).mark_line(
        clip=True,
        strokeWidth=2,
        opacity=1,
    ).encode(
        x=alt.X("dist:Q"),
        y=alt.Y('norm_density:Q'),
        tooltip=[alt.Tooltip("dist:Q", title="")]
    )
    combined_up = ecdf_plot_up_unselected+density_plot_up

    # Downstream
    ecdf_plot_down_unselected = alt.Chart(plot_df_down_0).mark_line(
        clip=True,
        strokeWidth=2,
        opacity=1,
    ).encode(
        x=alt.X("dist:Q",
                scale=alt.Scale(domain=[0, max_dist])),
        y=alt.Y('ecdf:Q', scale=alt.Scale(reverse=True, domain=[0,1])),
        color=alt.Color('snvscategorical:N', scale=None, 
        sort=alt.EncodingSortField('snvscategorical')),
        tooltip=[alt.Tooltip("dist:Q", title="")]
    )
    density_plot_down = alt.Chart(density_df_down).mark_line(
        clip=True,
        strokeWidth=2,
        opacity=1,
    ).encode(
        x=alt.X("dist:Q"),
        y=alt.Y('norm_density:Q'),
        tooltip=[alt.Tooltip("dist:Q", title="")]
    )
    combined_down = ecdf_plot_down_unselected+density_plot_down

    # Plot together
    combined = combined_up | combined_down
    combined.save(args.output_html)
    
  

if __name__=="__main__":
  main()