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
print(args)
entropy.df = read.csv(args$entropies, sep=' ', header=F, stringsAsFactors = F)
print(entropy.df)
if (args$relative=="F"){
  p = ggplot(entropy.df, aes(V1, V2))+
    geom_line()+
    theme_bw()+
    theme(panel.grid = element_blank())+
    ylab("Positional entropy (block diversity)")+
    xlab("Position (bp)")+
    theme(panel.border = element_blank())+
    theme(axis.line = element_line())
}
if (args$relative=="T"){
  head(entropy.df)
  p = ggplot(entropy.df, aes(V2, V3, colour=V1))+
    geom_line()+
    theme_bw()+
    theme(panel.grid = element_blank())+
    ylab("Positional entropy (block diversity)")+
    xlab("Distance from gene (bp)")+
    theme(panel.border = element_blank())+
    theme(axis.line = element_line())+
    scale_color_manual(values=c("red", "black"))+
    labs(colour="direction")
}

if (args$normalised=="Y"){
  p = p + ylim(c(0,1))
}


ggsave(p, file=args$output_pdf, width=as.numeric(args$width), height=as.numeric(args$height))
