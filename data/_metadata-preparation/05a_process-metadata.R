# Processes raw metadata downloaded from NCBI,
# - adding information on whether accession is chromosome/plasmid
# - formatting the dates properly

metadata <- read.csv('metadata-raw.csv',
                     header=T,
                     stringsAsFactors = F, 
                     row.names = 1,
                     na.strings="nan") # useful for reading in python-generated csv from pandas df
metadata$year <- as.numeric(gsub("-.*", "", metadata$collection_date))
chromosomes  <- read.csv('..//data-processing/CARD-chromosomes.csv', header=F, stringsAsFactors = F)$V1
plasmids  <- read.csv('../data-processing/CARD-plasmids.csv', header=F, stringsAsFactors = F)$V1

metadata$type <- sapply(rownames(metadata), 
                        function(x) ifelse(x %in% chromosomes, "chromosome", 
                                           ifelse(x %in% plasmids, "plasmid", "unknown")))

# Format dates
metadata$year = as.numeric(sapply(metadata$collection_date,
                                  function(x) gsub("-.*", "", x)))
metadata$public.year = as.numeric(sapply(metadata$first_public,
                                         function(x) gsub("-.*", "", x)))

metadata = metadata[order(metadata$biosample_id),]

write.csv(metadata, file=paste0("metadata-processed-v1.csv"))



