'Plot entropies. 

Usage:
  plot_entropies_cumulative.R <entropies> [--output_pdf=<pdf>] [--width=<w>] [--height=<h>] [--logscale=<l>]

Options:
  -h --help     Show this screen.
  --version     Show version.
  --output_pdf=<pdf> Name of output pdf. [default: plot_output.pdf]
  --width=<w> Width of output pdf. [default: 10]
  --height=<h> Height of output pdf. [default: 6]
  --logscale=<l> Whether to use logscale for y-axis. [default: Y]
# Entropies are assumed relative here!
' -> doc

suppressMessages(library(docopt, quietly = TRUE))
suppressMessages(library(ggplot2, quietly=TRUE))
suppressMessages(library(cowplot, quietly=TRUE))
suppressMessages(library(dplyr, quietly=TRUE))

args <- docopt(doc, version = 'Plot cumulative entropies v1.0')
entropy.df = read.csv(args$entropies, sep=' ', header=F, stringsAsFactors = F)
# Compute cumulative sum
entropy.df = entropy.df %>% group_by(V1, V2) %>% group_by(V1) %>% mutate(d=cumsum(V3))

  p = ggplot(entropy.df, aes(V2, d, colour=V1))+
    geom_line()+
    theme_bw()+
    theme(panel.grid = element_blank())+
    ylab("Cumulative entropy (block diversity)")+
    xlab("Distance from gene (bp)")+
    theme(panel.border = element_blank())+
    theme(axis.line = element_line())+
    scale_color_manual(values=c("red", "black"))+
    labs(colour="direction")

if (args$logscale=="Y"){
  p = p + scale_y_log10()
}



ggsave(p, file=args$output_pdf, width=as.numeric(args$width), height=as.numeric(args$height))
