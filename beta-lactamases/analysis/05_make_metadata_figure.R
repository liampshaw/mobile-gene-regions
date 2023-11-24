# Make plot of metadata
library(ggplot2)
library(tidyr)
library(dplyr)


theme_basic <- function () { 
  theme_bw(base_size=12) %+replace% 
    theme(
      axis.text=element_text(colour="black")
    ) %+replace% 
    theme(
      panel.grid=element_blank()
    )
}

GENUS_COUNT_THRESHOLD = 50
# Genus table
metadata_df = read.csv('../data/metadata.csv',
                       header=T,
                       stringsAsFactors = F,
                       row.names = 1)
genus.counts = sort(table(metadata_df$TaxGenus))
low.abundance.genera = names(genus.counts[which(genus.counts<GENUS_COUNT_THRESHOLD)])
# assign 'other' for these genera

GENERA = c("Other", names(genus.counts[which(genus.counts>GENUS_COUNT_THRESHOLD)]))
metadata_df$TaxGenus.simple = sapply(metadata_df$TaxGenus,
                                     function(x) ifelse(x %in% low.abundance.genera, "Other", x))
metadata_df$TaxGenus.simple = ordered(metadata_df$TaxGenus.simple, 
                                      levels=GENERA)
REGIONS = c("missing",
            "Sub-Saharan Africa", 
            "Middle East & North Africa",
            "Latin America & Caribbean",
            "South Asia",
            "East Asia & Pacific",
            "Europe & Central Asia",
            "North America")


genes = read.csv('../../genes.txt', header=F)$V1


supp_table = data.frame()
for (i in seq(1, length(genes))){
  gene = genes[i]
  gene_accs = read.csv(paste0('../data/accs-by-gene/', gene, '-accs.txt'), 
                              header=F)$V1
  metadata_df_gene = metadata_df[gene_accs,]
  d = metadata_df_gene %>% group_by(TaxGenus.simple, Year)%>%
    summarise(n=length(Year))
  d$Year = as.numeric(as.character(d$Year))
  plot_name =  paste0("p.", i)
   p = ggplot(d, aes(Year, n, fill=TaxGenus.simple))+
    geom_bar(stat="identity")+
    ggtitle(paste0(gene, " (n=", nrow(metadata_df_gene), ")"))+
    scale_fill_manual(values=c("grey", RColorBrewer::brewer.pal(n=11, name="Paired")),
                      drop=FALSE)+
    theme_basic()+
     theme(legend.position = "none")+
     ylab("")+
     xlab("")+
     xlim(c(2000,2022))
   assign(plot_name,p)
   contig_df = metadata_df_gene %>% group_by(Contig)%>%summarise(count=length(Contig))
   contig_df$Contig = ordered(contig_df$Contig,
                              levels=c("chromosome", "plasmid"))
   p.pie = ggplot(contig_df, aes(x = "", y = count, fill = Contig)) +
     geom_bar(width = 1, stat = "identity") +
     coord_polar(theta = "y") +
     theme_void()+
     theme(legend.position = "none")+
     scale_fill_manual(values=c("grey", "black"), drop=FALSE)
   assign(paste0("p.pie.", i), p.pie)
   # Summarise for supplementary table
   gene_table = metadata_df_gene %>% group_by(TaxGenus)%>%
     summarise(n.chromosome=length(Contig[Contig=="chromosome"]),
               n.plasmid=length(Contig[Contig=="plasmid"])) 
   gene_table$gene = gene
   gene_table = gene_table[,c("gene", "TaxGenus", "n.chromosome", "n.plasmid")]
   supp_table = rbind(supp_table, gene_table)
}

write.csv(supp_table, file='../../manuscript/figs/supp_table_contig_types.csv', 
          row.names = F,
          quote=F)

combined_plot = cowplot::plot_grid(p.1, p.2, p.3, p.4, 
                   p.5, p.6, p.7, p.8,
                   p.9, p.10, p.11, p.12, nrow=4, ncol=3, align='hv')
  
ggsave(cowplot::plot_grid(combined_plot), width=10,height=9,
       file="../../manuscript/figs/supp_fig_genus_counts.pdf")

combined_plot_pie = cowplot::plot_grid(p.pie.1, p.pie.2, p.pie.3, p.pie.4, 
                                       p.pie.5, p.pie.6, p.pie.7, p.pie.8,
                                       p.pie.9, p.pie.10, p.pie.11, p.pie.12, nrow=4, ncol=3)
ggsave(cowplot::plot_grid(combined_plot_pie), width=10,
       file="../../manuscript/figs/supp_fig_pie.pdf")

# legend = p.1 + theme(legend.position = "right")
# added in inkscape

# Do the same for chromosome/plasmid
for (i in seq(1, length(genes))){
  gene = genes[i]
  gene_accs = read.csv(paste0('../data/accs-by-gene/', gene, '-accs.txt'), 
                       header=F)$V1
  metadata_df_gene = metadata_df[gene_accs,]
  d = metadata_df_gene %>% group_by(TaxGenus.simple, Contig)%>%
    summarise(n=length(Year))
  plot_name =  paste0("p.contigs.", i)
  d$TaxGenus.simple = ordered(d$TaxGenus.simple,
                              levels=c(sort(GENERA[2:length(GENERA)]), "Other"))
  if (i>9){
    p = ggplot(d, aes(TaxGenus.simple,n, fill=Contig))+
      geom_bar(stat="identity", position=position_dodge(preserve="single"))+
      ggtitle(gene)+
      theme_basic()+
      scale_x_discrete(drop=FALSE)+
      theme(legend.position = "none")+
      ylab("")+
      scale_fill_manual(values=c("black", "red"))+
      theme(axis.text.x=  element_text(angle=45, hjust=1, face="italic")
)+
      xlab("")
    
  }
  else{
    p = ggplot(d, aes(TaxGenus.simple,n, fill=Contig))+
      geom_bar(stat="identity", position=position_dodge(preserve="single"))+
      ggtitle(gene)+
      theme_basic()+
      scale_x_discrete(drop=FALSE)+
      theme(legend.position = "none")+
      ylab("")+
      scale_fill_manual(values=c("black", "red"))+
      theme(axis.text.x=element_blank())+
      xlab("")
  }
  
  assign(plot_name,p)
}
combined_plot = cowplot::plot_grid(p.contigs.1, p.contigs.2, p.contigs.3, p.contigs.4, 
                                   p.contigs.5, p.contigs.6, p.contigs.7, p.contigs.8,
                                   p.contigs.9, p.contigs.10, p.contigs.11, p.contigs.12, nrow=4, ncol=3)


# Whether included in final 


# Genes by genus
# Summarise metadata with plots
library(ggplot2)
library(dplyr)

# Number of genes by genus
# >4 for no. of beta-lactamases
metadata_df$N.betalactamase.hits.simple = sapply(metadata_df$N.betalactamase.hits, 
                                                 function(x) ifelse(x>4, ">4", x))
# Summarise number of genes by genus
gene.hits.by.genus = metadata_df %>% group_by(TaxGenus.simple, N.betalactamase.hits.simple) %>%
  summarise(n=length(TaxGenus.simple))

# Sort by total number in genus
gene.hits.by.genus$TaxGenus.simple = ordered(gene.hits.by.genus$TaxGenus.simple,
                                             levels=c("Other", names(genus.counts[which(genus.counts>10)])))

# Summarise number of genes by genus and chromosome/plasmid
gene.hits.by.genus.chrom = metadata_df %>% group_by(TaxGenus.simple, Contig, N.betalactamase.hits.simple) %>%
  summarise(n=length(TaxGenus.simple))
gene.hits.by.genus.chrom$Contig[gene.hits.by.genus.chrom$Contig=="chromosome"] = "ncbi_chromosome"
gene.hits.by.genus.chrom$Contig[gene.hits.by.genus.chrom$Contig=="plasmid"] = "ncbi_plasmid"

# Order by total amounts
gene.hits.by.genus.chrom$TaxGenus.simple = ordered(gene.hits.by.genus.chrom$TaxGenus.simple,
                                                   levels=c("Other", names(genus.counts[which(genus.counts>GENUS_COUNT_THRESHOLD)])))
gene.hits.by.genus.chrom$N.betalactamase.hits.simple = ordered(gene.hits.by.genus.chrom$N.betalactamase.hits.simple,
                                                               levels=c("1", "2", "3", "4", ">4"))

gene.hits.by.genus.chrom$pseudo.n = gene.hits.by.genus.chrom$n+0.001



p.genes.genus = ggplot(gene.hits.by.genus.chrom, aes(TaxGenus.simple, n, fill=N.betalactamase.hits.simple))+
  geom_bar(stat="identity", position="stack", colour="black", size=0.1)+
  theme_basic()+
  coord_flip()+
  scale_y_continuous(breaks=seq(0,2500, 500),
                     minor_breaks=seq(0,2700,100),
                     position = "right")+
  theme(panel.grid.major.x = element_line(colour="grey"),
        panel.grid.minor.x = element_line(colour="grey"))+
  ylab("")+
  xlab("")+
  labs(fill="No. of\nbeta-lactamases\ncarried")+
  theme(axis.text.y=element_text(face="italic"))+
  scale_fill_brewer(palette="Reds")+
  theme(panel.border  = element_blank(),
        axis.line=element_line(colour="black"),
        axis.text.y=element_text(size=8))+
  facet_wrap(~Contig)+
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size=10))

# ggsave(file="../output/fig-genes-by-genus.pdf",
#        p.genes.genus,
#        width=8, height=4)
