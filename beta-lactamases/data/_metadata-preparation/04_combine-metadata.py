# Merge the dataframes back together to get a final 
# raw metadata file: metadata-raw.csv

import glob
import pandas as pd
import re

# These have been used in extract-xml.py - need to be the same
metadata_columns = ['organism', 'subspecies', 'strain', 'host', 'host_taxid', 'isolation_source', 'env_local_scale', 'collected_by',  'geo_loc_name', 'collection_date', 'first_public', 'lat_lon'] 
		

id_dict = {}
# Create dictionary of nuccore id and biosamples
with open('id_and_biosample.csv', 'r') as f:
	for i, line in enumerate(f.readlines()):
		nuccore_id, biosample_id = line.strip().split(',')
		id_dict[nuccore_id] = biosample_id

df_list = []
csv_files = glob.glob("tmp*xml.csv")
df_list = [None for x in range(len(csv_files))]
for i, csv_filename in enumerate(csv_files):
	print(csv_filename)
	df_list[i] = pd.read_csv(csv_filename, index_col=0)
# Combine the separate dataframes, dropping duplicates
df = pd.concat(df_list)
print(len(df))
unique_biosample_ids = df.index.drop_duplicates(keep ='first') 
print(len(unique_biosample_ids))

row_indexes = [list(df.index).index(x) for x in unique_biosample_ids]
df.index = range(len(df))


df = df.loc[row_indexes]
df.index = unique_biosample_ids
print(len(df))
#print(df.loc["SAMD00195990"])

# add the info to the id_dict
for nuccore_id in id_dict.keys():
	if id_dict[nuccore_id] in df.index: # if biosample info, add it
		biosample_id = id_dict[nuccore_id]
		metadata_list = [str(x) for x in list(df.loc[id_dict[nuccore_id]])]
		metadata_list.insert(0, biosample_id)
		id_dict[nuccore_id] = metadata_list	
	else:
		id_dict[nuccore_id] = ["" for x in range(len(metadata_columns)+1)] # no info available

#df.to_csv("metadata.csv")

#print(id_dict)

#print(id_dict["NZ_AP022513.1"])
print(id_dict["NZ_LT838197.1"])
metadata_columns.insert(0, "biosample_id")
with open("metadata-raw.csv", "w") as f:
	f.write("%s\n" % ",".join(metadata_columns))
	for nuccore_id, metadata in id_dict.items():
		#print(nuccore_id, metadata)
		f.write("\"%s\",%s\n" % (nuccore_id, re.sub(",nan", ",", ",".join(["\""+str(x)+"\"" for x in metadata]))))
