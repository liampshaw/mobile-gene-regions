# AUC summary
library(ggplot2)
theme_basic <- function () { 
  theme_bw(base_size=12) %+replace% 
    theme(
      axis.text=element_text(colour="black")
    ) %+replace% 
    theme(
      panel.grid=element_blank()
    )
}


auc_files = Sys.glob("AUCs/*deduplicated.txt")
summary_table = matrix(nrow=length(accession_files), ncol=6)
i = 1
all_results = matrix(nrow=0, ncol=4)
for (file in auc_files){
  gene = gsub("_.*", "", gsub("AUCs\\/", "", file))
  results = read.csv(file, header=F, sep=" ", stringsAsFactors = F)
  results = cbind(rep(gene, nrow(results)), results )
  all_results = rbind(all_results, results)
  #summary_table[i,] = c(gene, results)
  i = i + 1 
}
colnames(all_results) = c("gene", "snvs", "up", "down")
all_results = data.frame(all_results)
all_results$up = as.numeric(all_results$up)
all_results$down = as.numeric(all_results$down)

omit_genes = c("GES-24", "IMP-4", "OXA-10", "PER-1", "VIM-1")

all_results = all_results[which(!all_results$gene %in% omit_genes),]
all_results$snvs = ordered(all_results$snvs, levels=c(0, 1, 2, 3, 4, 5, 6, 7, ">7"))
all_results$snvs.ordinal = as.numeric(all_results$snvs)
p.up = ggplot(all_results, aes(snvs, up, fill=gene))+
  geom_line(aes(group=gene, colour=gene))+
  geom_point(shape=21, colour="black")+
  ylim(c(0,0.75))+
  ylab("AUC for pairwise breakpoints")+
  ggtitle("Upstream")+
  theme_basic()+
  theme(legend.position = "none")+
  scale_color_brewer(palette="Paired")+
  scale_fill_brewer(palette="Paired")+
  xlab("SNVs relative to focal gene")

p.down = ggplot(all_results, aes(snvs, down, fill=gene))+
  geom_line(aes(group=gene, colour=gene))+
  geom_point(shape=21, colour="black")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  coord_cartesian(ylim=c(0,0.75))+
  ylab("")+
  ggtitle("Downstream")+
  theme_basic()+
  scale_color_brewer(palette="Paired")+
  scale_fill_brewer(palette="Paired")+
  xlab("SNVs relative to focal gene")+
  labs(fill="Beta-lactamase")+
  guides(colour=FALSE)

p.combined = cowplot::plot_grid(p.up, p.down, rel_widths = c(1,1.5))
ggsave(p.combined, filename="../manuscript/figs/figure-breakpoint-AUCs.pdf",
       width=8, height=5)


# REPEAT FOR ALL 
