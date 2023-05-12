import altair as alt
import csv
import json
import pandas as pd
# from vega_datasets import data


def my_theme():
  return {
    'config': {
      'view': {'continuousHeight': 300, 'continuousWidth': 400},  # from the default theme
      'range': {'category': {'scheme': 'yelloworangered'}}
    }
  }
alt.themes.register('my_theme', my_theme)
alt.themes.enable('my_theme')


# Open the CSV file for reading
GENE = "CTX-M-15"
with open(GENE+'.output_dists.csv', newline='') as csvfile:
    # Read the CSV file into a dictionary
    reader = csv.DictReader(csvfile)

    # Convert the dictionary to a list of rows
    rows = list(reader)

# Convert the list of rows to a JSON string
json_str = json.dumps(rows)

# Print the JSON string
with open(GENE+'.output_dists.json', 'w') as f:
        f.write(json_str)

df = pd.read_json(GENE+'.output_dists.json')
print(df)
# Convert the data types of the columns to float if necessary
df['distup'] = df['dist.up'].astype(float)
df['distdown'] = df['dist.down'].astype(float)
# def map_to_categorical_snp(value):
#     if value <= 7:
#         return {0: chr(value), 1: chr(value), 2: chr(value), 3: chr(value), 4: chr(value), 5: chr(value), 6: chr(value), 7: chr(value)}
#     else:
#         return {">7": value}

def map_to_categorical_snp(values):
    return_values = [""] * len(values)
    for i,v in enumerate(values):
        if v<=7:
            return_values[i] = str(v)
        else:
            return_values[i] = ">7"
    return return_values

#print(list(df['snps']))
df['snpscategorical'] = map_to_categorical_snp(list(df['snps']))

max_dist = max(df['distup'])-10 # to clip the edges where it artifically goes to zero



cdf4_up = alt.Chart(df).transform_window(
    ecdf="cume_dist()",
    sort=[{"field": "distup"}],
).mark_line(
    interpolate="step-before", clip=True
).encode(
    x=alt.X("distup:Q",
        scale=alt.Scale(reverse=True, domain=[0,max_dist])),
    y=alt.Y('ecdf:Q', scale=alt.Scale(reverse=True)),
    color='snpscategorical:N'
)
cdf4_down = alt.Chart(df).transform_window(
    ecdf="cume_dist()",
    sort=[{"field": "distdown"}],
).mark_line(
    interpolate="step-before", clip=True
).encode(
    x=alt.X("distdown:Q", 
        scale=alt.Scale(domain=[0, max_dist])),
    y=alt.Y('ecdf:Q', scale=alt.Scale(reverse=True)),
    color='snpscategorical:N'
)
cdf4 = cdf4_up | cdf4_down
cdf4.save(GENE+'-cdf.html')


# Try adding dropdown selection
options = ['0', '1', '2', '3', '4', '5', '6', '7', '>7']
labels = [option + ' ' for option in options]

input_dropdown = alt.binding_radio(
    # Add the empty selection which shows all when clicked
    options=options + [None],
    labels=labels + ['All'],
    name='SNPs: '
)
selection = alt.selection_point(
    fields=['snpscategorical'],
    bind=input_dropdown
)
# input_dropdown = alt.binding_select(options=["0","1","2","3","4","5","6","7",">7"], 
#     name='SNPs ')
# selection = alt.selection_point(fields=['snpscategorical'], bind=input_dropdown)
color = alt.condition(
    selection,
    alt.Color('snpscategorical:N').legend(None),
    alt.value('lightgray')
)

cdf5_up = alt.Chart(df).transform_window(
    ecdf="cume_dist()",
    sort=[{"field": "distup"}],
).mark_line(
    interpolate="step-before", clip=True
).encode(
    x=alt.X("distup:Q",
        scale=alt.Scale(reverse=True, domain=[0,max_dist])),
    y=alt.Y('ecdf:Q', scale=alt.Scale(reverse=True)),
    color=color
)#.add_params(
  #  selection
#).transform_filter(
#    selection
#)



cdf5_down = alt.Chart(df).transform_window(
    ecdf="cume_dist()",
    sort=[{"field": "distdown"}],
).mark_line(
    interpolate="step-before", clip=True
).encode(
    x=alt.X("distdown:Q", 
        scale=alt.Scale(domain=[0, max_dist])),
    y=alt.Y('ecdf:Q', scale=alt.Scale(reverse=True)),
    #color='snpscategorical:N'
    color=color
)#.add_params(
 #   selection
#).transform_filter(
#    selection
#)



legend = alt.Chart(df).mark_point().encode(
    alt.Y('snpscategorical:N').axis(orient='right'),
    color=color
).add_params(
    selection
)

cdf5_up_selected = alt.Chart(df).transform_filter(
    selection
).transform_window(
    ecdf="cume_dist()",
    sort=[{"field": "distup"}],
).mark_line(
    interpolate="step-before",
    clip=True,
    strokeWidth=2,
    opacity=1,
).encode(
    x=alt.X("distup:Q",
            scale=alt.Scale(reverse=True, domain=[0, max_dist])),
    y=alt.Y('ecdf:Q', scale=alt.Scale(reverse=True, domain=[0,1])),
    color=alt.Color('snpscategorical:N').legend(None),
)

cdf5_up_unselected = alt.Chart(df).transform_filter(
    ~selection
).transform_window(
    ecdf="cume_dist()",
    sort=[{"field": "distup"}],
).mark_line(
    interpolate="step-before",
    clip=True,
    strokeWidth=2,
    opacity=0.4,
).encode(
    x=alt.X("distup:Q",
            scale=alt.Scale(reverse=True, domain=[0, max_dist])),
    y=alt.Y('ecdf:Q', scale=alt.Scale(reverse=True, domain=[0, 1])),
     color=alt.Color('snpscategorical:N', scale=None#, 
#        sort=alt.EncodingSortField('snpscategorical')),
)
)

cdf5_down_selected = alt.Chart(df).transform_filter(
    selection
).transform_window(
    ecdf="cume_dist()",
    sort=[{"field": "distdown"}],
).mark_line(
    interpolate="step-before",
    clip=True,
    strokeWidth=2,
    opacity=1,
).encode(
    x=alt.X("distdown:Q",
            scale=alt.Scale(domain=[0, max_dist])),
    y=alt.Y('ecdf:Q', scale=alt.Scale(reverse=True, domain=[0, 1])),
    color=alt.Color('snpscategorical:N').legend(None),
)

cdf5_down_unselected = alt.Chart(df).transform_filter(
    ~selection
).transform_window(
    ecdf="cume_dist()",
    sort=[{"field": "distdown"}],
).mark_line(
    interpolate="step-before",
    clip=True,
    strokeWidth=2,
    opacity=0.4,
).encode(
    x=alt.X("distdown:Q",
            scale=alt.Scale(domain=[0, max_dist])),
    y=alt.Y('ecdf:Q', scale=alt.Scale(reverse=True, domain=[0,1])),
     color=alt.Color('snpscategorical:N', scale=None#, 
#        sort=alt.EncodingSortField('snpscategorical')),
)
)

cdf5_up = (cdf5_up_unselected+cdf5_up_selected ).add_selection(selection)
cdf5_down = (cdf5_down_unselected+cdf5_down_selected ).add_selection(selection)



# alt.layer(
#   cdf5_up.add_params(selection),
#   cdf5_up.transform_filter(selection).encode( f"toDate(datum.snpscategorical) == {selection.snpscategorical}")
# )

cdf5_combined = cdf5_up | cdf5_down | legend
cdf5_combined.save(GENE+'-cdf-selector.html')

cdf5_down.save(GENE+'-cdf-selector-down.html')



