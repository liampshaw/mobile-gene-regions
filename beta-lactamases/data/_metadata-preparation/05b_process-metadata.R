# Reading in processing and edited metadata
# and format further 
# for final saving

library(countrycode) # for countries
library(AMR) # for taxonomy

# For processing metadata after manual editing. 

metadata.edited = read.csv("metadata-processed-v2.csv",
                           row.names = 1,
                           header=T,
                           stringsAsFactors = F)

# Add the number of gene hits
metadata.edited$N.betalactamase.hits = sapply(metadata.edited$gene.hits, function(x) length(strsplit(x, ",")[[1]]))



metadata.edited$TaxSpecies = sapply(metadata.edited$organism, 
                                    function(x) paste( strsplit( x, split =" ")[[1]][1], 
                                                       strsplit( x, split =" ")[[1]][2]))

metadata.edited$TaxGenus = gsub(" .*", "", metadata.edited$TaxSpecies)






metadata.edited$TaxFamily = AMR::mo_family(metadata.edited$TaxSpecies)
# Assign Stutzerimonas as Pseudomonadaceae (see https://www.sciencedirect.com/science/article/pii/S0723202021001120?via%3Dihub)
metadata.edited$TaxFamily[which(metadata.edited$TaxSpecies=="Stutzerimonas stutzeri")] = "Pseudomonadaceae"

# Add order
metadata.edited$TaxOrder = AMR::mo_order(metadata.edited$TaxFamily)
# Add class
metadata.edited$TaxClass = AMR::mo_class(metadata.edited$TaxFamily)
# Add phylum
metadata.edited$TaxPhylum = AMR::mo_phylum(metadata.edited$TaxFamily)

################
# COUNTRY
################
metadata.edited$Country = metadata.edited$country.name
metadata.edited$Country_ISO = countrycode(metadata.edited$Country, 
                                          origin="country.name",
                                          destination = "iso3c")

# Assign region
metadata.edited$Region = countrycode(metadata.edited$Country_ISO,
                                     origin="iso3c",
                                     destination="region")
# If NA, missing
metadata.edited$Region[is.na(metadata.edited$Region)] = "missing"


metadata.edited$Year = metadata.edited$year
metadata.edited$BioSampleID = metadata.edited$biosample_id
metadata.edited$NucCoreID = rownames(metadata.edited)


metadata.edited$NCBI.host = metadata.edited$host
metadata.edited$NCBI.isolation_source = metadata.edited$isolation_source
metadata.edited$NCBI.env_local_scale = metadata.edited$env_local_scale
metadata.edited$NCBI.collected_by = metadata.edited$collected_by

metadata.edited$Gene.hits = metadata.edited$gene.hits

metadata.edited$Contig = metadata.edited$type 

# Column names I want
use.names = c("NucCoreID",
              "BioSampleID", 
              "N.betalactamase.hits",
              "Gene.hits",
              "Country",
              "Region",
              "Year",
              "TaxSpecies",
              "TaxGenus", 
              "TaxFamily", 
              "TaxOrder", 
              "TaxClass", 
              "TaxPhylum",
              "NCBI.host", 
              "NCBI.isolation_source",
              "NCBI.env_local_scale",
              "NCBI.collected_by",
              "Contig")


# Write everything to file
write.csv(metadata.edited[,use.names], file="metadata-processed-v3.csv")

