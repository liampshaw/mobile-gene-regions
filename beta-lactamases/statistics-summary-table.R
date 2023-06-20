# Create statistics on the beta-lactamase datasets

metadata = read.csv("../data/metadata.csv", 
                    header=T, 
                    stringsAsFactors = F,
                    row.names = 1)

summariseMetadata <- function(accs_file){
  accs = read.csv(accs_file, header=F, stringsAsFactors = F)$V1
  metadata.reduced = metadata[accs,]
  earliest.genome = min(na.omit(as.numeric(metadata.reduced$Year)))
  number.total = nrow(metadata.reduced)
  number.chromosome = length(metadata.reduced$Contig[which(metadata.reduced$Contig=="chromosome")])
  number.plasmid = length(metadata.reduced$Contig[which(metadata.reduced$Contig=="plasmid")])
  number.genera = length(table(metadata.reduced$TaxGenus))
  
  return(c(earliest.genome, number.total, number.chromosome, number.plasmid, number.genera))
}

accession_files <- sort(Sys.glob("accs/*.txt"))
summary_table = matrix(nrow=length(accession_files), ncol=6)
i = 1
for (file in accession_files){
  gene = gsub("_.*", "", gsub("accs\\/", "", file))
  results = summariseMetadata(file)
  summary_table[i,] = c(gene, results)
  i = i + 1 
}
summary.table = data.frame(summary_table)
colnames(summary.table) = c("Gene", "Earliest genome", "N (total)", "N (chrom)", "N (plas)", "N (genera)")
write.csv(file="summary_table.csv", summary.table, quote=F, row.names = F)

