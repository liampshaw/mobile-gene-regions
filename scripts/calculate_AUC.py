import csv
import json
import pandas as pd
import argparse

def get_options():
    parser = argparse.ArgumentParser(description='Calculate AUC of the empirical cdf of breakpoint distances')
    parser.add_argument('--dist_csv', help='CSV of output distances, produced from pangraph pipeline.', type=str)
    parser.add_argument('--gene_of_interest', help='Gene of interest to subset to comparisons involving only this named variant', type=str, default='', required=False)   
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

def auc(x, y, max_dist=5000):
    '''returns the area under breakpoint distance curve as specified by a set of x, y points
    N.B. we do 1-auc because of defining it as the empirical cdf which goes to 1'''
    x = list(x)
    y = list(y)
    x.append(max_dist)
    y.append(1)
    diffs_x = [x[i+1]-x[i] for i in range(len(x)-1)]
    return(1-sum([diffs_x[i] * y[i+1] for i in range(1, len(diffs_x))])/max_dist)

def main(): 
    # Open the CSV file for reading
    args = get_options()

    if args.gene_of_interest!='':
        gene = args.gene_of_interest.upper()
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

     # Convert the data types of the columns to float 
    df['distup'] = df['dist.up'].astype(float)
    df['distdown'] = df['dist.down'].astype(float)
    # And add categorical SNVs
    df['snvscategorical'] = map_to_categorical_snv(list(df['snvs']))
    plot_df_up = df.groupby('snvscategorical').apply(lambda x: ecdf(x['distup'])).reset_index(level=1, drop=True)
    plot_df_up = plot_df_up.assign(snvscategorical=plot_df_up.index)
    plot_df_down = df.groupby('snvscategorical').apply(lambda x: ecdf(x['distdown'])).reset_index(level=1, drop=True)
    plot_df_down = plot_df_down.assign(snvscategorical=plot_df_down.index)

    # summarise by SNV level
    for x in sorted(set(plot_df_up["snvscategorical"])):
        print(x, 
            auc(plot_df_up["dist"].loc[plot_df_up["snvscategorical"]==x], 
        plot_df_up["ecdf"].loc[plot_df_up["snvscategorical"]==x]),
            auc(plot_df_down["dist"].loc[plot_df_down["snvscategorical"]==x], 
        plot_df_down["ecdf"].loc[plot_df_down["snvscategorical"]==x]))


if __name__=="__main__":
    main()

