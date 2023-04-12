'Plot pangenome blocks linearly with annotations from a gff.

Usage:
  plot_with_gff_annotations.R <block_csv> <gff_file> <focal_gene> <focal_block> [--output_pdf=<pdf>] [--width=<w>] [--height=<h>] 

Options:
  -h --help     Show this screen.
  --version     Show version.
  --output_pdf=<pdf> Name of output pdf. [default: plot_output.pdf]
  --width=<w> Width of output pdf. [default: 8]
  --height=<h> Height of output pdf. [default: 4]
  N.B. currently assumes that annotations are for only a single genome (hence the genome block requirement)
' -> doc

require(ggplot2, quietly=TRUE)
require(cowplot, quietly=TRUE)
require(gggenes, quietly=TRUE)
require(ggdendro, quietly=TRUE)
require(reshape2, quietly=TRUE)
suppressPackageStartupMessages(library(vegan, quietly = TRUE))
require(docopt, quietly=TRUE)
options(warn=-1)



args <- docopt(doc, version = 'Plot with gff annotations v1.0')
GENE = args$focal_gene
GFF.FILE = args$gff_file

genome.blocks = read.csv(args$block_csv, header=T, stringsAsFactors = F)
genome.blocks$forward <- ifelse(genome.blocks$strand=="+", TRUE, FALSE)
block.counts <- table(genome.blocks$block)
blocks.which.need.colours <- names(block.counts)[which(block.counts>1)]

genome.blocks$block.coloured <- sapply(genome.blocks$block,
                                       function(x) ifelse(x %in% blocks.which.need.colours,
                                                          x,
                                                          "_other"))
block.colours <- unique(genome.blocks$colour)
names(block.colours) <- unique(genome.blocks$block.coloured)


FOCAL.BLOCK = read.csv(args$focal_block, header=F, stringsAsFactors=F)$V1
gene.block.locations <- genome.blocks[which(genome.blocks$block==FOCAL.BLOCK),c("start", "end")]
rownames(gene.block.locations) <- genome.blocks[which(genome.blocks$block==FOCAL.BLOCK),"genome"]
transformAnnotationsBlocks <- function(genome, default_offset=9926){
  offset <- default_offset-gene.block.locations[genome, c("start")]
  return(offset)
}

genome.blocks$new.start <- apply(genome.blocks, MARGIN=1, function(x) as.numeric(x["start"])+transformAnnotationsBlocks(x["genome"]))-10001
genome.blocks$new.end <- apply(genome.blocks, MARGIN=1, function(x) as.numeric(x["end"])+transformAnnotationsBlocks(x["genome"]))-10001

# Get all the genome paths (in terms of blocks)
genome.paths <- sapply(unique(genome.blocks$genome), function(x) paste(genome.blocks[which(genome.blocks$genome==x), "block"], collapse=","))
genome.paths <- sort(genome.paths, decreasing = TRUE)

# Subset to unique paths
# Keep one representative genome for each
genome.path.reps <- names(genome.paths)[!duplicated(genome.paths)]
#genome.blocks.unique <- genome.blocks[which(genome.blocks$genome %in% genome.path.reps),]
# And store the number of examples of each
#genome.blocks.unique$genome.path <- genome.paths[genome.blocks.unique$genome]
#genome.blocks.unique$n.reps <- sapply(genome.blocks.unique$genome.path,
 #                                     function(x) table(genome.paths)[x])
#genome.blocks.unique$genome.path.name <- paste0("Type", as.numeric(as.factor(genome.blocks.unique$genome.path )))
#genome.blocks.unique$genome.n <- sapply(genome.blocks.unique$n.reps, function(x) ifelse(x==1,
                #                                                                        "", paste0("n=", x)))

# Now use a dendrogram to order the genomes
m <- acast(genome ~ block, data=genome.blocks, fill=0, fun.aggregate=length)
m.dist <- vegdist(m, method="jaccard") # jaccard distances based on block presence/absence
dendro <- as.dendrogram(hclust(m.dist))
dendro_order <- order.dendrogram(dendro)
genome_labels <- dendro_data(dendro)$labels$label
genome.blocks$genome.ordered <- ordered(genome.blocks$genome,
                                               levels=genome_labels)


annotation.hits = read.csv(GFF.FILE, 
                        sep='\t', 
                        header=F, 
                        comment.char = "#", 
                        stringsAsFactors = F)
annotation.hits = annotation.hits[which(annotation.hits$V3=="CDS"),]
annotation.hits$product = gsub(" ", "\n", gsub(";.*", "", gsub(".*product=", "", annotation.hits$V9)))

GENOMES = unique(annotation.hits$V1)

# Convert to positive strand for the focal gene
convertAnnotations <- function(input_df, focal_gene){
  combined_df = NULL
  for (genome in unique(input_df$V1)){
    df = input_df[which(input_df$V1==genome),]
    focal.product = df$product[grep(focal_gene, df$product)]
    if (length(focal.product)!=0){
focal.gene.strand = df[which(df$product==focal.product),"V7"]
    if (focal.gene.strand=="+"){
      focal.gene.start = df$V4[grep(focal_gene, df$product)]
      focal.gene.end = df$V5[grep(focal_gene, df$product)]
      df$new.start = df$V4-focal.gene.start
      df$new.end = df$V5-focal.gene.start
      df$new.strand = df$V7
    }
    if (focal.gene.strand=="-"){
      # We are going to reverse-complement everything
      # Choose an arbitrary maximum length D
      # to reverse things from
      # 0...10---20........500 where gene on -ve strand 10-20
      # becomes after reverse-complementing
      # 0.................480+++490...500
      # i.e. we use a maximum length to offset and flip strands
      D = max(df$V5)
      df$new.start = D-df$V5
      df$new.end = D-df$V4
      df$new.strand = ifelse(df$V7=="-", "+", "-")
      focal.gene.start = df$new.start[grep(focal_gene, df$product)]
      df$new.start = df$new.start-focal.gene.start
      df$new.end = df$new.end-focal.gene.start
    }
    df$forward = ifelse(df$new.strand=="+", T, F)
    if (is.null(combined_df)){
      combined_df = df
    }
    else{
      combined_df = rbind(combined_df, df)
    }
    }
    
  }
  return(combined_df)
}
# 

annotation.hits = convertAnnotations(annotation.hits, GENE)

# Reduce to only those with annotations provided
#annotation.hits = annotation.hits[which(annotation.hits$V1 %in% GENOMES),]

#genome.blocks = genome.blocks[which(genome.blocks$genome.ordered %in% GENOMES),]
#genome.blocks$genome.ordered =ordered(genome.blocks$genome.ordered, )

limits = c(min(c(genome.blocks$new.start,genome.blocks$new.end)),
       max(c(genome.blocks$new.start,genome.blocks$new.end)))

annotation.hits$genome.ordered = ordered(annotation.hits$V1,
                                         levels=levels(genome.blocks$genome.ordered))
annotation.hits = annotation.hits[which(annotation.hits$new.start>limits[1]*1.1 & annotation.hits$new.end<limits[2]*1.1),]

# Genome plot

p.genome = ggplot(genome.blocks[which(genome.blocks$genome.ordered %in% GENOMES),], 
                  aes(xmin = new.start, xmax = new.end, forward = forward, 
                    y = genome.ordered, fill = block.coloured)) +
  geom_gene_arrow(arrow_body_height = unit(6, "mm"),
                  arrowhead_height = unit(0, "mm"),
                  arrowhead_width = unit(0, "mm"),
                  alpha=1)+
  theme_genes()+
  scale_fill_manual(values=block.colours)+
  ylab("")+
  theme(legend.position = "none")+
  scale_y_discrete(breaks=levels(genome.blocks$genome.ordered), labels=levels(genome.blocks$genome.ordered))+
  theme(plot.title=element_text(hjust=0.5))+
  xlab("Position (bp)")

genome_levels = levels(genome.blocks$genome.ordered)
genome_levels_reduced = genome_levels[which(genome_levels %in% GENOMES)]
bottom.genome = genome_levels_reduced[1]
top.genome = genome_levels_reduced[length(genome_levels_reduced)]

p.genome = p.genome+
      expand_limits(y= c(0, length(levels(as.numeric(genome.blocks$genome.ordered))) + 2))

p.final = p.genome+
  ggrepel::geom_text_repel(aes(label=block,
            x=(new.start+new.end)/2, 
            y=genome.ordered), nudge_y = -0.25,
            size=1,
            data=genome.blocks[which(genome.blocks$genome.ordered==bottom.genome),])+
  geom_gene_arrow(aes( xmin=new.start, xmax=new.end, 
                                        y=genome.ordered,
                                        forward=forward), 
                                   inherit.aes = FALSE, 
                                   size=0.5,
                                   annotation.hits,
                                   arrow_body_height = unit(2, "mm"),
                                   arrowhead_height = unit(2, "mm"),
                                   arrowhead_width = unit(1, "mm"),
                                   alpha=1) +
  ggrepel::geom_text_repel(aes(label=product, 
                               x=(new.start+new.end)/2, 
                               y=genome.ordered),
                           size=1,
                           data=annotation.hits[which(annotation.hits$genome.ordered==top.genome),], inherit.aes = FALSE,
                           hjust=0,
                           min.segment.length =0, force = 10, nudge_y = 0.5)


pdf(args$output_pdf, width=as.numeric(args$width), height=as.numeric(args$height))
p.final
dev.off()

