'Plot a tree of beta-lactamase variants in a dataset. 

Usage:
  plot_tree_variants.R --aln=<aln> --tree=<tree> --variants=<variants> --dup_names=<dup_names> [--output_pdf=<pdf>] [--width=<w>] [--height=<h>] 

Options:
  -h --help          Show this screen.
  --version          Show version.
  <aln>              Alignment file of the variants (DNA).
  <tree>             Newick tree file to be plotted (produced from alignment file).
  <variants>         Assignments of isolates to named protein variants.
  <dup_names>        File with names of duplicated sequences (produced from seqkit rmdup).
  --output_pdf=<pdf> Name of output pdf. [default: plot_output.pdf]
  --width=<w>        Width of output pdf. [default: 8]
  --height=<h>       Height of output pdf. [default: 8]

' -> doc

library(docopt, quietly = TRUE)
library(ape, quietly=TRUE)
library(ggplot2, quietly=TRUE)
suppressMessages(library(ggtree, quietly=TRUE))
library(RColorBrewer, quietly=TRUE)
options(warn=-1)

args <- docopt(doc, version = 'Plot entropies v1.0')


variant_alignment_file = args$aln
tree_file = args$tree
tree = read.tree(tree_file)
root=NULL


dup.names = read.csv(args$dup_names, header=F, stringsAsFactors=F, sep='\t')
dup.names$rep = gsub(",.*", "", dup.names$V2)
rownames(dup.names) = dup.names$rep



variants = read.csv(args$variants, header=F, stringsAsFactors = F, row.names = 1)

if (!is.null(root)){
  tree = root(tree, outgroup=root)
}

dna <- read.dna(variant_alignment_file, format='fasta')
isolates = gsub(" .*", "",labels(dna))
snvs = as.numeric(sapply(labels(dna), function(x) strsplit(x, " ", fixed=TRUE)[[1]][4]))



variant.names = variants[isolates, "V2"]
n.values = dup.names[isolates,"V1"]
n.values[is.na(n.values)] = 1 # if they only have one rep, not in dedup.txt file, so must be added
  metadata = data.frame(isolate=isolates, 
                      variant=variant.names,
                      n=n.values,
                      snvs=snvs)

rownames(metadata) <- metadata$isolate

gene.family = unique(gsub("-.*$", "", variant.names[grep("-", variant.names)]))

metadata$variant[which(metadata$variant=="unnamed")] = paste0(gene.family, "-?")
metadata$variant[which(metadata$variant=="truncated")] = "*"

metadata$isolate2 = metadata$isolate
metadata$n = as.numeric(as.character(metadata$n))

metadata$tip.label = paste0(metadata$variant, " (n=", metadata$n, ")")

p <- ggtree(tree, layout = 'rectangular')

# Add colors
snvs.categorical.colour.palette <- RColorBrewer::brewer.pal(name="RdYlBu", n=9)
metadata$snvs.categorical = sapply(metadata$snvs, function(x) ifelse(x<8, x, ">7"))
metadata$snvs.categorical = ordered(metadata$snvs.categorical, levels=c(seq(0,7), ">7"))

p <- p %<+% metadata+
  geom_tippoint(aes(size=n, fill=snvs.categorical), shape=21, colour="black")+
  geom_tiplab(aes(label=tip.label), size=3, offset = 0.2)+
  theme(legend.title= element_text(size=8), 
        legend.text = element_text(size=6))+
  scale_size_continuous(range=c(2, 10), breaks=c(1, 10, 50, 100), name="No. sequences")+
  geom_treescale(width = 1, offset = 0.1, linesize = 1)+
  guides(fill="none")+
  scale_fill_manual(values=snvs.categorical.colour.palette)
expanded.range = as.numeric(layer_scales(p)$x$range$range[2]*1.1)
p = p +xlim(c(0, expanded.range)) # expand plot so as to not obscure labels



pdf(args$output_pdf, width=8, height=8)
p
dev.off()

