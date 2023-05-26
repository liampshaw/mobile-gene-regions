'Plot entropies. 

Usage:
  plot_entropies.R <entropies> [--output_pdf=<pdf>] [--width=<w>] [--height=<h>] [--normalised=<n>] [--relative=<r>]

Options:
  -h --help     Show this screen.
  --version     Show version.
  --output_pdf=<pdf> Name of output pdf. [default: plot_output.pdf]
  --width=<w> Width of output pdf. [default: 10]
  --height=<h> Height of output pdf. [default: 6]
  --normalised=<n> Whether entropies are normalised to 1. [default: Y]
  --relative=<r> Whether distances are upstream/downstream relative to focal gene. [default: F]

' -> doc

suppressMessages(library(docopt, quietly = TRUE))
suppressMessages(library(ggplot2, quietly=TRUE))
suppressMessages(library(cowplot, quietly=TRUE))
suppressMessages(library(dplyr, quietly=TRUE))

args <- docopt(doc, version = 'Plot entropies v1.0')
entropy.df = read.csv(args$entropies, sep=' ', header=F, stringsAsFactors = F)
if (args$relative=="F"){
  p = ggplot(entropy.df, aes(V1, 1-V2))+
    geom_line()+
    theme_bw()+
    theme(panel.grid = element_blank())+
    ylab("Positional entropy (block diversity)")+
    xlab("Position (bp)")+
    theme(panel.border = element_blank())+
    theme(axis.line = element_line())
}
if (args$relative=="T"){
  dist.max <- max(entropy.df$V2)
  p.together = ggplot(entropy.df, aes(V2, 1-V3, colour=V1))+
    geom_line()+
    theme_bw()+
    theme(panel.grid = element_blank())+
    ylab("Positional entropy (block diversity)")+
    xlab("Distance from gene (bp)")+
    theme(panel.border = element_blank())+
    theme(axis.line = element_line())+
    scale_color_manual(values=c("red", "black"))+
    labs(colour="direction")
    p.upstream = ggplot(entropy.df[which(entropy.df$V1=="upstream"),], aes(-V2, 1-V3))+
    geom_line()+
    theme_bw()+
    theme(panel.grid = element_blank())+
    ylab("1 - positional entropy (block diversity)")+
    xlab("Distance from gene (bp)")+
    theme(panel.border = element_blank())+
    theme(axis.line = element_line())+
    scale_y_continuous(limits=c(0,1), expand = c(0, 0))+
    scale_x_continuous(limits=c(-dist.max, 0), expand=c(0,0))+
    annotate(geom="segment", x=0, y=0,xend=0, yend=1,  linetype='dashed')
    p.downstream = ggplot(entropy.df[which(entropy.df$V1=="downstream"),], aes(V2, 1-V3))+
    geom_line()+
    theme_bw()+
    theme(panel.grid = element_blank())+
    ylab("Positional entropy (block diversity)")+
    xlab("Distance from gene (bp)")+
    theme(panel.border = element_blank())+
    theme(axis.line = element_line())+
    scale_y_continuous(limits=c(0,1), expand = c(0, 0))+
    scale_x_continuous(limits=c(0, dist.max), expand=c(0,0))+    
    theme(axis.line.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          axis.title.y=element_blank())+
    annotate(geom="segment", x=0, y=0,xend=0, yend=1, linetype='dashed')
    p = cowplot::plot_grid(p.upstream, p.downstream, nrow=1)
}

if (args$normalised=="Y"){
  p = p + ylim(c(0,1))
}


ggsave(p, file=args$output_pdf, width=as.numeric(args$width), height=as.numeric(args$height))
