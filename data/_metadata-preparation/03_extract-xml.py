# Extract the metadata as csv from the downloaded xml files

import untangle
import sys
import pandas as pd
import re

xml_file = sys.argv[1]
#xml_file = 'combined.xml'
# obj = untangle.parse(xml_file)
#
# biosample_attributes = list(obj.DocumentSummarySet.DocumentSummary.SampleData.BioSample.Attributes.Attribute)
# attribute_dict = {a.get_attribute('attribute_name'): a.cdata for a in biosample_attributes}
# organism = obj.DocumentSummarySet.DocumentSummary.Organism.cdata
# attribute_dict['organism'] = organism
# print(attribute_dict)

# For reading a whole xml file (new method)
with open(xml_file, 'r', encoding='utf-8') as f:
    newItemBool = True
    xml_list = []
    xml_index = 0
    for i, line in enumerate(f.readlines()):
        if i>3: # skip header...
            #print(line)
            if line.strip()=="</DocumentSummarySet>":
                print('end')
                xml_index = 'end'
            if line.strip()=="<DocumentSummary>":
                xml_list.append('')
            if xml_index!='end':
                xml_list[xml_index] += line.encode("utf-8", "ignore").decode()
            if line.strip()=="</DocumentSummary>":
                xml_index += 1


biosample_dict = {}
biosample_dict_harmonized = {}
for i, xml_string in enumerate(xml_list):
    print(i)
    xml_string_ascii = xml_string.encode("ascii", "ignore").decode("ascii", "ignore") # to ignore unicode characters e.g. accents in submitter names
    try:
        #obj = untangle.parse(str(xml_string.encode("utf-8", "ignore")))
        obj = untangle.parse(xml_string_ascii)
    except Exception:
        print(xml_string_ascii)
        pass #print('error, skipping')
    biosample_accession = re.sub('BioSample:', '', str(obj.DocumentSummary.SourceSample.cdata))
    #print(biosample_accession)
    biosample_attributes = list(obj.DocumentSummary.SampleData.BioSample.Attributes.Attribute)
    #print(biosample_attributes)
    harmonized_attribute_dict = {a.get_attribute('harmonized_name'): a.cdata for a in biosample_attributes} # only take harmonized names
    attribute_dict = {a.get_attribute('attribute_name'): a.cdata for a in biosample_attributes} # save unharmonized names just in case
    # add in organism info
    organism = obj.DocumentSummary.Organism.cdata
    harmonized_attribute_dict['organism'] = organism
    attribute_dict['organism'] = organism
    # add in 'first public' dates, if they exist
    first_public = {a.get_attribute('attribute_name'): a.cdata for a in biosample_attributes if 'public' in a.get_attribute('attribute_name').lower()}
    #type_material = {a.get_attribute('type-material'): a.cdata for a 
    # Get the earliest of these dates
    first_public_earliest = next(iter(sorted(first_public.values())), None)
    harmonized_attribute_dict['first_public'] = first_public_earliest
    #print(first_public_earliest)
    biosample_dict_harmonized[biosample_accession] = harmonized_attribute_dict
    biosample_dict[biosample_accession] = attribute_dict
    #print(i, biosample_accession, attribute_dict)

#
print(biosample_dict_harmonized)
# This is my current list of useful metadata
useful_columns = ['organism', 'subspecies', 'strain', 'host', 'host_taxid', 'isolation_source', 'env_local_scale', 'collected_by',  'geo_loc_name', 'collection_date', 'first_public', 'lat_lon'] 
#useful_columns = ['organism', 'host', 'collected by', 'geographic location', 'collection date', 'latitude and longitude']
# 'ENA first public', 'INSDC first public'

df = pd.DataFrame.from_dict(biosample_dict_harmonized, orient='index')

# Add empty columns if needed
empty_columns = [x for x in useful_columns if x not in df.columns]
df2 = df.reindex(columns = df.columns.tolist() + empty_columns)
df2[useful_columns].to_csv(xml_file+'.csv')

# Also save unharmonized
full_df = pd.DataFrame.from_dict(biosample_dict, orient='index')
empty_columns = [x for x in useful_columns if x not in full_df.columns]
full_df2 = df.reindex(columns = full_df.columns.tolist() + empty_columns)
full_df2.to_csv(xml_file+'.full.csv') # also save all the metadata
# create object for each of the separated files - read from <DocumentSummary> to </DocumentSummary>

