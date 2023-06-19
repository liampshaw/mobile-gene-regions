'Plot output distances from pangraph output.

Usage:
  plot_output_dists.R <dists> <dedup_names> <variant_assignments> [--output_pdf_prefix=<prefix>] [--focal_gene=<g>]

Options:
  -h --help     Show this screen.
  --version     Show version.
  --output_pdf_prefix=<prefix> Prefix for output pdfs [default: plot_output]
  --focal_gene=<g> Name of focal gene variant

' -> doc

require(docopt, quietly = TRUE)
require(ggplot2, quietly=TRUE)
require(ggridges, quietly=TRUE)
require(cowplot, quietly=TRUE)
require(RColorBrewer, quietly=TRUE)
require(ape, quietly=TRUE)
suppressMessages(require(ggtree, quietly=TRUE))
require(gridExtra, quietly=TRUE)
#require(magick)

args <- docopt(doc, version = 'Plot_output_dists v1.0')




d <- read.csv(args$dists,
              header=T,
              stringsAsFactors = F)

d$snps.categorical <- sapply(d$snps, function(x)
  ifelse(x<8, x, ">7"))
d$snps.categorical <- ordered(d$snps.categorical,
                      levels=c(seq(0,7), ">7"))
snps.categorical.colour.palette <- RColorBrewer::brewer.pal(name="RdYlBu", n=9)



# Also get most common gene sequence - assuming it exists
metadata <- read.csv(args$dedup_names, 
                     sep='\t', 
                     header=F, 
                     col.names = c("n", "isolates"))
most.abundant.variant = strsplit(metadata$isolates[1], split=", ")[[1]]

d.subset.compare.to.most.common = d[which(d$seq1 %in% most.abundant.variant | d$seq2 %in% most.abundant.variant),]

# Compare to focal gene sequence
gene.assignments = read.csv(args$variant_assignments, header=F, stringsAsFactors=F)
variants = gene.assignments$V2
names(variants) = gene.assignments$V1
focal.gene.isolates = gene.assignments$V1[which(gene.assignments$V2==args$focal_gene)]
d.subset.compare.to.focal = d[which(d$seq1 %in% focal.gene.isolates | d$seq2 %in% focal.gene.isolates),]



dist.max <- max(c(d$dist.up, d$dist.down))


makePlots <- function(df){
    # Add number of comparisons for each SNV distance
    snp.comparison.table = data.frame(SNVs=names(table(df$snps.categorical)),
                                    N.pairs=as.numeric(table(df$snps.categorical)))
    table.theme <- gridExtra::ttheme_default(
    core = list(fg_params=list(cex = 0.5)),
    colhead = list(fg_params=list(cex = 0.5)),
    rowhead = list(fg_params=list(cex = 0.5)))

    p.upstream <-  ggplot(df, aes( -dist.up, group=snps.categorical, colour=snps.categorical))+
      stat_ecdf()+
      theme_bw()+
    annotate(geom="segment", x=0, y=0,xend=0, yend=1,  linetype='dashed')+
      theme(legend.position = "none")+
      labs(colour="SNVs")+
      ylab("cdf")+
      xlab("distance from gene (bp)")+
      theme(panel.grid = element_blank())+
      scale_color_manual(values=snps.categorical.colour.palette)+
      theme(axis.line.x = element_line(colour = "black"),
        axis.line.y=element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        title = element_text())+
      scale_y_continuous(limits=c(0,1), expand = c(0, 0))+
      scale_x_continuous(limits=c(-dist.max, 0), expand=c(0,0))+
      theme(plot.title=element_text(hjust=0.5))+
      annotation_custom(tableGrob(snp.comparison.table, theme=table.theme, rows=NULL), 
                        xmax=-dist.max*3.5/5, ymax=1.5)


    p.downstream <- ggplot(df, aes( dist.down,group=snps.categorical, colour=snps.categorical))+
      stat_ecdf()+
      theme_bw()+
    annotate(geom="segment", x=0, y=0,xend=0, yend=1,  linetype='dashed')+
      labs(colour="SNVs")+
      ylab("cdf")+
      xlab("distance from gene (bp)")+
      theme(panel.grid = element_blank())+
      scale_color_manual(values=snps.categorical.colour.palette)+
      scale_y_continuous(limits=c(0,1), expand = c(0, 0))+
      scale_x_continuous(limits=c(0, dist.max), expand=c(0,0))+      
      theme(axis.line.y=element_blank(),
            axis.text.y=element_blank(),
            axis.title.y=element_blank())+
      ggtitle("")+
      theme(axis.line.x = element_line(colour = "black"),
        axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
        theme(plot.title=element_text(hjust=0.5))+
      scale_y_reverse()

    p = cowplot::plot_grid(p.upstream, p.downstream, 
                           rel_widths  = c(0.855, 1),
                           align = "h")


    return(p)
}


pdf(paste0(args$output_pdf_prefix, '-all.pdf'), width=8, height=4)
makePlots(d)
dev.off()

pdf(paste0(args$output_pdf_prefix, '-compare-to-most-common.pdf'), width=8, height=4)
makePlots(d.subset.compare.to.most.common)
dev.off()

pdf(paste0(args$output_pdf_prefix, '-compare-to-focal-gene.pdf'), width=8, height=4)
makePlots(d.subset.compare.to.focal)
dev.off()



