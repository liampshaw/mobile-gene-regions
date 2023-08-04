regular.packages = c("ape", "cowplot", "docopt", "ggdendro", "gridExtra",
             "gggenes", "reshape2", "RColorBrewer", "vegan")

for (p in regular.packages){
  if (!require(p, quietly = TRUE))
    install.packages(p, repos = "http://cran.us.r-project.org")
}

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")

