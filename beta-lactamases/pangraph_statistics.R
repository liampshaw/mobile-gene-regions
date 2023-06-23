stats = read.csv('gfastats/combined_stats.csv', header=T, stringsAsFactors = F)

stats$blocks.per.seq = stats$blocks/stats$sequences
stats$circular_segments.per.seq = stats$circular_segments/stats$sequences
stats$bubbles.per.seq = stats$bubbles/stats$sequences
#stats$total_size.per.seq = stats$total_size/stats$sequences
stats$fasta_size.per.seq = stats$total_fasta_size/stats$sequences


library(tidyr)

stats.per.seq = stats[,grep(".per.seq", colnames(stats))]
stats.per.seq$gene = stats$gene
stats.per.seq$fasta_size = stats$total_fasta_size
stats.per.seq$average_block_length = stats$average_block_length
stats.per.seq$block_N50 = stats$block_N50

colnames(stats.per.seq) = c("N_blocks (per seq)", "Circular segments (per seq)", 
                            "Bubbles (per seq)", "Pangraph size (per seq)", 
                            "gene", "fasta_size", "Average block length (bp)",
                            "Block length N50")

stats.df = stats.per.seq %>% pivot_longer(cols=c("N_blocks (per seq)", 
                                                 "Circular segments (per seq)", 
                                                 "Bubbles (per seq)", 
                                                 "Pangraph size (per seq)",
                                                 "Average block length (bp)",
                                                 "Block length N50"))
ggplot(stats.df, aes(gene, value))+
  geom_bar(stat="identity", fill="black")+
  facet_wrap(~name, scale="free")+
  theme_basic()+
  ylab("")+
  xlab("")+
  coord_flip()

