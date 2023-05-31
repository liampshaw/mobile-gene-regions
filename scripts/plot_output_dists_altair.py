import altair as alt
import csv
import json
import pandas as pd
import argparse

def my_theme():
  return {
    'config': {
      'view': {'continuousHeight': 300, 'continuousWidth': 400},  # from the default theme
      'range': {'category': {'scheme': 'yelloworangered'}}
    }
  }
alt.themes.register('my_theme', my_theme)
alt.themes.enable('my_theme')




interpolate_option = "step"

def get_options():
    parser = argparse.ArgumentParser(description='Produce html plot of pairwise breakpoint distances')
    parser.add_argument('--dist_csv', help='CSV of output distances, produced from pangraph pipeline.', type=str)
    parser.add_argument('--output_html', help='Output html (default: output.html)', type=str, default='output.html')
    parser.add_argument('--gene_of_interest', help='Gene of interest to subset to comparisons involving only this named variant', type=str, default='', required=False)   
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
  
def map_to_categorical_snp(values):
    '''Maps SNPs to categories - have >7 as (arbitrary) highest category'''
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

def main(): 
    # Open the CSV file for reading
    args = get_options()

    if args.gene_of_interest!='':
        df = pd.read_csv(args.dist_csv)
        df = df.loc[(df["variant1"]==args.gene_of_interest) | (df["variant2"]==args.gene_of_interest)]
        tmp_csv = args.dist_csv+'.'+args.gene_of_interest+'.csv'
        df.write_csv(tmp_csv)
        # Convert csv to json
        tmp_json = tmp_csv+'_tmp.json'
        convert_csv_to_json(csv_input=tmp_csv, json_output=tmp_json)
        df = pd.read_json(tmp_json)

    else:
        # Convert csv to json
        tmp_json = args.dist_csv+'_tmp.json'
        convert_csv_to_json(csv_input=args.dist_csv, json_output=tmp_json)
        df = pd.read_json(tmp_json)

    # Convert the data types of the columns to float 
    df['distup'] = df['dist.up'].astype(float)
    df['distdown'] = df['dist.down'].astype(float)
    # And add categorical SNPs
    df['snpscategorical'] = map_to_categorical_snp(list(df['snps']))
    # Maximum distance
    # to clip the edges where ecdf artifically goes to 1 (so plot drops to 0) because we don't look beyond
    max_dist = max(df['distup'])-10 # 

    # For dropdown selection for SNP category of comparison
    options = sorted(set(df['snpscategorical']))
    labels = [option + ' ' for option in options]
    # Sort labels in the same order
    labels_sorted = [label for _, label in sorted(zip(options, labels))]

    input_dropdown = alt.binding_radio(
        options=options + [None],
        labels=labels_sorted + ['All'],
        name='SNPs: '
    )

    selection = alt.selection_point(
        fields=['snpscategorical'],
        bind=input_dropdown,
        #empty='all'
    )

    color = alt.condition(
        selection ,# & (alt.datum.snpscategorical != 'None'),
        alt.Color('snpscategorical:N').legend(None),
        alt.value('lightgray')
    )



    # Create ecdf dataframes
    plot_df_up = df.groupby('snpscategorical').apply(lambda x: ecdf(x['distup'])).reset_index(level=1, drop=True)
    plot_df_up = plot_df_up.assign(snpscategorical=plot_df_up.index)
    plot_df_down = df.groupby('snpscategorical').apply(lambda x: ecdf(x['distdown'])).reset_index(level=1, drop=True)
    plot_df_down = plot_df_down.assign(snpscategorical=plot_df_down.index)


    # The strategy for selection is that 
    # we use an unselected one as a background layer 
    # and put the selected one on top

    # Produce legend
    legend = alt.Chart(plot_df_up).mark_point().encode(
        alt.Y('snpscategorical:N').axis(orient='right'),
        color=color
    ).add_params(
        selection
    )

    # Upstream plots
    ecdf_plot_up_selected = alt.Chart(plot_df_up).transform_filter(
        selection
    ).mark_line(
        interpolate=interpolate_option,
        clip=True,
        strokeWidth=2,
        opacity=1,
    ).encode(
        x=alt.X("dist:Q",
                scale=alt.Scale(reverse=True, domain=[0, max_dist])),
        y=alt.Y('ecdf:Q', scale=alt.Scale(reverse=True, domain=[0,1])),
        color=alt.Color('snpscategorical:N').legend(None),
    )
    ecdf_plot_up_unselected = alt.Chart(plot_df_up).mark_line(
        interpolate=interpolate_option,
        clip=True,
        strokeWidth=2,
        opacity=0.8,
    ).encode(
        x=alt.X("dist:Q",
                scale=alt.Scale(reverse=True, domain=[0, max_dist])),
        y=alt.Y('ecdf:Q', scale=alt.Scale(reverse=True, domain=[0,1])),
        color=alt.Color('snpscategorical:N', scale=None, 
        sort=alt.EncodingSortField('snpscategorical')),
        tooltip=[alt.Tooltip("dist:Q", title="")]
    )


    # Downstream plots
    ecdf_plot_down_selected = alt.Chart(plot_df_down).transform_filter(
        selection
    ).mark_line(
        interpolate=interpolate_option,
        clip=True,
        strokeWidth=2,
        opacity=1,
    ).encode(
        x=alt.X("dist:Q",
                scale=alt.Scale(domain=[0, max_dist])),
        y=alt.Y('ecdf:Q', scale=alt.Scale(reverse=True, domain=[0, 1])),
        color=alt.Color('snpscategorical:N').legend(None),
    )
    ecdf_plot_down_unselected = alt.Chart(plot_df_down).mark_line(
        interpolate=interpolate_option,
        clip=True,
        strokeWidth=2,
        opacity=0.8,
    ).encode(
        x=alt.X("dist:Q",
                scale=alt.Scale(domain=[0, max_dist])),
        y=alt.Y('ecdf:Q', scale=alt.Scale(reverse=True, domain=[0,1])),
        color=alt.Color('snpscategorical:N', scale=None, 
        sort=alt.EncodingSortField('snpscategorical')),
        tooltip=[alt.Tooltip("dist:Q", title="")]
    )

    # Combine the plots
    ecdf_plot_up = (ecdf_plot_up_unselected+ecdf_plot_up_selected ).add_params(selection)
    ecdf_plot_down = ecdf_plot_down_unselected.add_params(selection)+ecdf_plot_down_selected 
    ecdf_plot_combined = ecdf_plot_up | ecdf_plot_down | legend

    # Save them
    ecdf_plot_combined.save(args.output_html)

if __name__=="__main__":
    main()
