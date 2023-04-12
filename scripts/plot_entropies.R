'Plot entropies. 

Usage:
  plot_entropies.R <entropies> [--output_pdf=<pdf>] [--width=<w>] [--height=<h>] 

Options:
  -h --help     Show this screen.
  --version     Show version.
  --output_pdf=<pdf> Name of output pdf. [default: plot_output.pdf]
  --width=<w> Width of output pdf. [default: 8]
  --height=<h> Height of output pdf. [default: 8]

' -> doc

library(docopt, quietly = TRUE)
library(ggplot2)
library(cowplot)
library(dplyr)

args <- docopt(doc, version = 'Plot entropies v1.0')

entropy.df = read.csv(args$entropies, sep=' ', header=F, stringsAsFactors = F)
p = ggplot(entropy.df, aes(V1, V2))+
    geom_line()+
    theme_bw()+
    theme(panel.grid = element_blank())+
    ylab("Positional entropy (block diversity)")+
    xlab("Position")+
    theme(panel.border = element_blank())+
    theme(axis.line = element_line())

ggsave(p, file=args$output_pdf, width=as.numeric(args$width), height=as.numeric(args$height))
