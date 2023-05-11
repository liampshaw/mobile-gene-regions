import altair as alt
import csv
import json
import pandas as pd
# from vega_datasets import data

# source = data.movies.url
# print(source)

# chart = alt.Chart(source).transform_window(
#     cumulative_count="count()",
#     sort=[{"field": "IMDB_Rating"}],
# ).mark_area().encode(
#     x="IMDB_Rating:Q",
#     y="cumulative_count:Q"
# )

# chart.save('chart.html')


# source = data.cars()

# brush = alt.selection_interval()
# points = alt.Chart(source).mark_point().encode(
#     x='Horsepower',
#     y='Miles_per_Gallon',
#     color=alt.condition(brush, 'Origin', alt.value('lightgray'))
# ).add_params(
#     brush
# )

# bars = alt.Chart(source).mark_bar().encode(
#     y='Origin',
#     color='Origin',
#     x='count(Origin)'
# ).transform_filter(
#     brush
# )

# b = points & bars
# b.save('both.html')


# import csv
# import json
# import pandas as pd

# # Open the CSV file for reading
# with open('GES-24.output_dists.csv', newline='') as csvfile:
#     # Read the CSV file into a dictionary
#     reader = csv.DictReader(csvfile)

#     # Convert the dictionary to a list of rows
#     rows = list(reader)

# # Convert the list of rows to a JSON string
# json_str = json.dumps(rows)

# # Print the JSON string
# with open('GES-24.output_dists.json', 'w') as f:
#         f.write(json_str)

# df = pd.read_json('GES-24.output_dists.json')

# # Convert the data types of the columns to float if necessary
# df['distup'] = df['dist.up'].astype(float)
# df['distdown'] = df['dist.down'].astype(float)

# # Create the Altair chart
# brush = alt.selection_interval()

# scatter = alt.Chart(df).mark_point().encode(
#     x='distup',
#     y='distdown',
#  color=alt.condition(brush, 'snps', alt.value('lightgray'))
# ).add_params(
#     brush
# )

# cdf = alt.Chart(df).transform_window(
#     cumulative_count='count()',
#     sort=[{"field": "distup"}],
# ).mark_area().encode(
#     x="distup:Q",
#     y="cumulative_count:Q",
#     color='snps'
# )

# print(df)
# d = pd.crosstab(df.distup, columns=df.snps).cumsum()
# d = d.stack().reset_index()
# print(d)
# d = d.rename(columns={0:'CummulativeCount'})
# d['snps'] = d['snps'].astype(object)

# print(d)
# cdf = alt.Chart(
#     d,
#     width=120,
#     height=80
# ).transform_window(
#     cumulative_count='count()',
#     sort=[{"field": "distup"}],
# ).transform_density(
#     'cumulative_count',
#     groupby=['snps'],
#     as_=['distup', 'density']
# ).mark_area().encode(
#     x="distup:Q",
#     y='density:Q',
#     tooltip='snps:N'
# )
# cdf2 = alt.Chart(d).mark_line().encode(x='distup:T', y='CummulativeCount:Q', color='snps')


# chart.save('GES-24-scatter.html')
# cdf2.save('GES-24-cdf.html')

# import numpy as np
# from numpy import *
# sq=df['distup'].value_counts()
# print(sq.sort_index().cumsum()*1./len(sq))

# # Compute the ECDF
# x = np.sort(df['distup'])
# y = np.arange(0, len(df)) / len(df)
# data_up = pd.DataFrame({'x':x, 'y':1-y})

# x_down = np.sort(df['distdown'])
# y_down = np.arange(0, len(df)) / len(df)
# data_down = pd.DataFrame({'x':x_down, 'y':1-y_down})

# #plt.plot(x, y, marker='.')
# # plt.show()
# cdf_up = alt.Chart(data_up
#     ).mark_line(
#     ).encode(x=alt.X('x:Q', scale=alt.Scale(reverse=True)), 
#     y='y:Q')
# cdf_down = alt.Chart(data_down).mark_line().encode(x='x:Q', y='y:Q')

# cdf_combined = cdf_up | cdf_down
# cdf_combined.save('GES-24-cdf3.html')



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
    bind=input_dropdown,
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
).add_params(
    selection
).transform_filter(
    selection
)

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
    color=alt.Color('snpscategorical:N').scale(domain=options),
).add_params(
    selection
).transform_filter(
    selection
)

cdf5_combined = cdf5_up | cdf5_down
cdf5_combined.save(GENE+'-cdf-selector.html')

cdf5_down.save(GENE+'-cdf-selector-down.html')


