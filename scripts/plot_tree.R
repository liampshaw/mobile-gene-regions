library(ape, quietly=TRUE)
suppressMessages(library(ggtree, quietly=TRUE))
library(ggplot2, quietly=TRUE)

FUNCTIONAL.CATEGORIES = c("Betalactamase", "Broad-spectrum","Carbapenemase","ESBL")
FUNCTIONAL.COLOURS = c("#a6cee3","#1f78b4","#b2df8a","#33a02c")
names(FUNCTIONAL.COLOURS) = FUNCTIONAL.CATEGORIES

args = commandArgs(trailingOnly=TRUE)

metadata.all = read.csv("data/BLDB-subset.tsv", sep="\t", header=T)
metadata.all$Functional.information[which(metadata.all$Functional.information=="")] = "Betalactamase"
#tree = read.tree(snakemake@input[[1]])
tree = read.tree(args[1])
p = ggtree(tree)+theme_tree2()

expanded.range = as.numeric(layer_scales(p)$x$range$range[2]*1.1)
print(expanded.range)
p = p +xlim(c(0, expanded.range)) # expand plot so as to not obscure labels
print(expanded.range)

rownames(metadata.all) = metadata.all$Protein.name
metadata = metadata.all[tree$tip.label,]
metadata$Isolate = metadata$Protein.name
metadata = metadata[,c("Isolate", "Functional.information")]



#pdf(snakemake@output[[1]], height=nrow(metadata)/10, width=8)
pdf(args[2], height=nrow(metadata)/10, width=8)
p %<+% metadata +geom_tiplab(aes(colour=Functional.information), size=2)+
	geom_tippoint(aes(colour=Functional.information))+
	scale_color_manual(values=FUNCTIONAL.COLOURS)
dev.off()
