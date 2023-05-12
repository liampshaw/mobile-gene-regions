import altair as alt
import csv
import json
import pandas as pd

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


# Open the CSV file for reading
GENE = "CTX-M-15"
# BELOW ONLY NEEDED FOR NOT YET CONVERTED CSV TO JSON
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
# ^^^ ONLY NEEDED FOR NOT YET CONVERTED CSV TO JSON


df = pd.read_json(GENE+'.output_dists.json')
# Convert the data types of the columns to float if necessary
df['distup'] = df['dist.up'].astype(float)
df['distdown'] = df['dist.down'].astype(float)

# Map SNPs to categories - have >7 as (arbitrary) highest category
def map_to_categorical_snp(values):
    return_values = [""] * len(values)
    for i,v in enumerate(values):
        if v<=7:
            return_values[i] = str(v)
        else:
            return_values[i] = ">7"
    return return_values

df['snpscategorical'] = map_to_categorical_snp(list(df['snps']))


max_dist = max(df['distup'])-10 # to clip the edges where it artifically goes to zero


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

# Produce legend
legend = alt.Chart(df).mark_point().encode(
    alt.Y('snpscategorical:N').axis(orient='right'),
    color=color
).add_params(
    selection
)


def ecdf(values):
    '''returns the ecdf for values'''
    v_sorted = [0]+sorted(values)
    return pd.DataFrame({'ecdf':[x/len(values) for x in range(0, len(values)+1)],
                    'dist':v_sorted})

# Create ecdf dataframes
plot_df_up = df.groupby('snpscategorical').apply(lambda x: ecdf(x['distup'])).reset_index(level=1, drop=True)
plot_df_up = plot_df_up.assign(snpscategorical=plot_df_up.index)
plot_df_down = df.groupby('snpscategorical').apply(lambda x: ecdf(x['distdown'])).reset_index(level=1, drop=True)
plot_df_down = plot_df_down.assign(snpscategorical=plot_df_down.index)

# The strategy for selection is that 
# we use an unselected one as a background layer 
# and put the selected one on top

# Upstream plots
cdf_plot_up_selected = alt.Chart(plot_df_up).transform_filter(
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
cdf_plot_up_unselected = alt.Chart(plot_df_up).mark_line(
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
cdf_plot_down_selected = alt.Chart(plot_df_down).transform_filter(
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
cdf_plot_down_unselected = alt.Chart(plot_df_down).mark_line(
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
cdf_plot_up = (cdf_plot_up_unselected+cdf_plot_up_selected ).add_params(selection)
cdf_plot_down = cdf_plot_down_unselected.add_params(selection)+cdf_plot_down_selected 
cdf_plot_combined = cdf_plot_up | cdf_plot_down | legend

# Save them
cdf_plot_combined.save(GENE+'-cdf-selector.html')
