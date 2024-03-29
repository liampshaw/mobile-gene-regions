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


auc_files = Sys.glob("AUCs/*AUCs_focal_gene_deduplicated.txt")
summary_table = matrix(nrow=length(auc_files), ncol=6)
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
all_results$snvs = ordered(all_results$snvs, levels=c(0, 1, 2, 3, 4, 5, 6, 7, ">7"))
all_results$snvs.ordinal = as.numeric(all_results$snvs)

all_results$auc.total = (all_results$up + all_results$down)/2 * 10000

# unsure whether to include this - just to make the figure nicer
# in terms of colours
omit_genes = c("GES-24", "IMP-4", "OXA-10", "PER-1", "VIM-1")
all_results_top = all_results[which(!all_results$gene %in% omit_genes),]

p.auc.overall = ggplot(all_results_top, aes(snvs, auc.total, colour=gene, fill=gene))+
  ylim(c(0,5000))+
  ylab("Average flanking region overlap (bp)")+
  theme_basic()+
  #geom_point(shape=21,colour="black")+
  scale_color_brewer(palette="Paired")+
  scale_fill_brewer(palette="Paired")+
  xlab("SNVs relative to focal gene")+
  labs(fill="Beta-lactamase")+
  guides(colour=FALSE)+
  theme(legend.position = "none")+
  stat_smooth(aes(x=snvs.ordinal-1, auc.total), method="lm", se=FALSE, fullrange = FALSE)+
  scale_x_continuous(breaks=seq(0,8), labels=c("0", "1", "2", "3", "4", "5", "6", "7", ">7"))
p.auc.facet = ggplot(all_results_top, aes(snvs, auc.total, colour=gene, fill=gene))+
  geom_path(aes(group=gene))+
  geom_point(shape=21, aes(fill=gene), colour="black")+
  ylab("AUC for breakpoint distance distribution")+
  theme_basic()+
  scale_color_brewer(palette="Paired")+
  scale_fill_brewer(palette="Paired")+
  xlab("SNVs relative to focal gene")+
  labs(fill="Beta-lactamase")+
  guides(colour=FALSE)+
  facet_wrap(~gene, scales="free_y")+
  theme(legend.position = "none")+
  theme(axis.text=element_text(size=12),
        strip.text = element_text(size=14))+
  ylab("")+xlab("")
p.auc.combined = cowplot::plot_grid(p.auc.overall+ggtitle("(a)")+theme(title=element_text(size=16)), 
                                    p.auc.facet+ggtitle("(b)")+theme(title=element_text(size=16)), align='h', 
                                    rel_widths = c(1,1.8),
                                    axis='v')

ggsave(p.auc.combined, 
       file='../../manuscript/figs/FINAL-figure-AUC-top-bl-final.pdf',
       width=9.5, height=5)
# 
cor.test(all_results$auc.total, as.numeric(all_results$snvs), method="spearman")

all_results_top = all_results[which(!all_results$gene %in% omit_genes),]
p.up = ggplot(all_results_top, aes(snvs, up, fill=gene))+
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

p.down = ggplot(all_results_top, aes(snvs, down, fill=gene))+
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

p.combined = cowplot::plot_grid(p.up, p.down, rel_widths = c(1,2))
ggsave(p.combined, filename="../../manuscript/figs/FINAL-figure-breakpoint-AUCs-top-genes.pdf",
       width=8, height=5)

p.up = ggplot(all_results, aes(snvs, up))+
  geom_line(aes(group=gene))+
  geom_point(shape=21, colour="black", fill="white")+
  ylim(c(0,0.75))+
  ylab("AUC for pairwise breakpoints")+
  ggtitle("Upstream")+
  theme_basic()+
  theme(legend.position = "none")+
  xlab("SNVs relative to focal gene")

p.down = ggplot(all_results, aes(snvs, down))+
  geom_line(aes(group=gene))+
  geom_point(shape=21, colour="black",fill="white")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  coord_cartesian(ylim=c(0,0.75))+
  ylab("")+
  ggtitle("Downstream")+
  theme_basic()+
  xlab("SNVs relative to focal gene")+
  theme(legend.position = "none")

p.combined = cowplot::plot_grid(p.up, p.down, rel_widths = c(1,1.1))
ggsave(p.combined, filename="../../manuscript/figs/FINAL-figure-breakpoint-AUCs-all.pdf",
       width=8, height=5)


# may need to REPEAT FOR ALL (i.e. not deduplicated) for supplementary

# genes cluster below the line of equality i.e.
# greater structural diversity upstream?
p.auc.upstream.downstream =ggplot(all_results[all_results$snvs==0,], aes(down, up))+
  geom_abline(intercept = 0, slope=1, linetype='dashed')+
  geom_point()+
  coord_fixed()+
  ylim(c(0,1))+
  xlim(c(0,1))+
  ggrepel::geom_text_repel(aes(label=gene))+
  theme_basic()+
  xlab("Downstream AUC")+
  ylab("Upstream AUC")
ggsave(p.auc.upstream.downstream, filename="../../manuscript/figs/FINAL-figure-AUC-upstream-downstream.pdf",
       width=6, height=6)
# 
wilcox.test(all_results[all_results$snvs==0,"down"]-all_results[all_results$snvs==0,"up"])


# Just non-integron-associated
all_results = all_results[which(!all_results$gene %in% c("GES-1", "VIM-1", "OXA-10", "IMP-4")),]
wilcox.test(all_results[all_results$snvs==0,"down"]-all_results[all_results$snvs==0,"up"])
