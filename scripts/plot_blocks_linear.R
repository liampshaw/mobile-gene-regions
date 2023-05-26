'Plot pangenome blocks in a linear way.

Usage:
  plot_blocks_linear.R <block_csv> [--focal_block=<bl>] [--output_pdf=<pdf>] [--width=<w>] [--height=<h>] [--strains=<s>]

Options:
  -h --help     Show this screen.
  --version     Show version.
  --focal_block=<bl>  Name of focal block file. [default: none].
  --output_pdf=<pdf> Name of output pdf. [default: plot_output.pdf]
  --width=<w> Width of output pdf. [default: 8]
  --height=<h> Height of output pdf. [default: 8]
  --strains=<s> File of strains to include. [default: all]

' -> doc

library(docopt, quietly = TRUE)
require(ggplot2, quietly = TRUE)
require(cowplot, quietly = TRUE)
require(gggenes, quietly = TRUE)
require(ggdendro, quietly = TRUE)
require(reshape2, quietly = TRUE)
suppressPackageStartupMessages(library(vegan, quietly = TRUE))
options(warn=-1)



args <- docopt(doc, version = 'Plot-blocks-linear v1.0')


# Genome blocks approach
genome.blocks <- read.csv(args$block_csv, header=T, stringsAsFactors = F)
#colnames(genome.blocks) <- c("genome", "block", "strand", "start", "end", "colour")
genome.blocks$forward <- ifelse(genome.blocks$strand=="+", TRUE, FALSE)
block.counts <- table(genome.blocks$block)
blocks.which.need.colours <- names(block.counts)[which(block.counts>1)]

genome.blocks$block.coloured <- sapply(genome.blocks$block,
                                       function(x) ifelse(x %in% blocks.which.need.colours,
                                                          x,
                                                          "_other"))
block.colours <- unique(genome.blocks$colour)
names(block.colours) <- unique(genome.blocks$block.coloured)

# If focal block
if (args$focal_block!='none'){
  focal_block = read.csv(args$focal_block)$V1[1]
}
if (args$focal_block=='none'){
  # Find maximum length
  max.length = max(genome.blocks$end)
  # Pick first possible block at centre
  central.point = floor(max.length/2)
  # Distance from centre for blocks in random genome
  first.genome = genome.blocks[which(genome.blocks$genome==unique(genome.blocks$genome)[1]),]
  dist_from_centre = abs(first.genome$start-central.point)
  focal_block = first.genome$block[which(dist_from_centre==min(dist_from_centre))]
}


#gene.block.locations <- genome.blocks[which(genome.blocks$block==focal_block),c("start", "end")]
#rownames(gene.block.locations) <- genome.blocks[which(genome.blocks$block==focal_block),"genome"]
#transformAnnotationsBlocks <- function(genome, default_offset=9926){
#  offset <- default_offset-gene.block.locations[genome, c("start")]
#  return(offset)
#}

focal_block.locations <- genome.blocks[which(genome.blocks$block==focal_block),c("start", "end")]
rownames(focal_block.locations) <- genome.blocks[which(genome.blocks$block==focal_block),"genome"]

# Need to add a different anchor block location for genomes that lack the anchor block
genomes.without.focal_block = as.character(unique(genome.blocks$genome)[which(!unique(genome.blocks$genome) %in% rownames(focal_block.locations))])
for (g in genomes.without.focal_block){
  focal_block.locations[g,] = c(0,0)
}

transformAnnotationsBlocks <- function(genome, default_offset=0){
  offset <- (default_offset+focal_block.locations[genome, c("start")]) 
  return(offset)
}

genome.blocks$new.start <- apply(genome.blocks, MARGIN=1, function(x) as.numeric(x["start"])-transformAnnotationsBlocks(x["genome"]))
genome.blocks$new.end <- apply(genome.blocks, MARGIN=1, function(x) as.numeric(x["end"])-transformAnnotationsBlocks(x["genome"]))

# Subset to the given strains
if (args$strains!="all"){
  strains = read.csv(args$strains, header=F, stringsAsFactors = F)$V1
  if (length(strains[which(strains %in% unique(genome.blocks$genome))])<2){
    stop("Your strains file must contain at least two strains in the pangraph.")
  }
  genome.blocks = genome.blocks[which(genome.blocks$genome %in% strains),]
}



# Get all the genome paths (in terms of blocks)
genome.paths <- sapply(unique(genome.blocks$genome), function(x) paste(genome.blocks[which(genome.blocks$genome==x), "block"], collapse=","))
genome.paths <- sort(genome.paths, decreasing = TRUE)

# Subset to unique paths in the given strains
# Keep one representative genome for each
genome.path.reps <- names(genome.paths)[!duplicated(genome.paths)]
genome.blocks.unique <- genome.blocks[which(genome.blocks$genome %in% genome.path.reps),]
# And store the number of examples of each
genome.blocks.unique$genome.path <- genome.paths[genome.blocks.unique$genome]
genome.blocks.unique$n.reps <- sapply(genome.blocks.unique$genome.path,
                                      function(x) table(genome.paths)[x])
genome.blocks.unique$genome.path.name <- paste0("Type", as.numeric(as.factor(genome.blocks.unique$genome.path )))
genome.blocks.unique$genome.n <- sapply(genome.blocks.unique$n.reps, function(x) ifelse(x==1,
                                                                                        "", paste0("n=", x)))

# Now use a dendrogram to order the genomes
m <- acast(genome ~ block, data=genome.blocks.unique, fill=0, fun.aggregate=length, value.var="genome.n")
m.dist <- vegdist(m, method="jaccard") # jaccard distances based on block presence/absence
if (length(m.dist)>0){ # if there is more than one block pattern!
  dendro <- as.dendrogram(hclust(m.dist))
  dendro_order <- order.dendrogram(dendro)
  genome_labels <- dendro_data(dendro)$labels$label
}
if (length(m.dist)==0){
  genome_labels=genome.blocks.unique$genome[1]
}

genome.blocks.unique$genome.ordered <- ordered(genome.blocks.unique$genome,
                                               levels=genome_labels)

# Make the linear block plot
p.blocks <- ggplot(genome.blocks.unique, aes(xmin = new.start, xmax = new.end, forward = forward, y = genome.ordered, fill = block.coloured)) +
  geom_gene_arrow(arrow_body_height = unit(2, "mm"),
                  arrowhead_height = unit(2, "mm"),
                  arrowhead_width = unit(1, "mm")) +
  theme_genes()+
  scale_fill_manual(values=block.colours)+
  ylab("")+
  theme(legend.position = "none")+
  scale_y_discrete(breaks=genome.blocks.unique$genome.ordered, labels=genome.blocks.unique$genome.n)+
  theme(plot.title=element_text(hjust=0.5))+
  xlab("Position (bp)")


p.final = p.blocks


pdf(file=args$output_pdf, width=as.numeric(args$width), height=as.numeric(args$height))
p.final
dev.off()

