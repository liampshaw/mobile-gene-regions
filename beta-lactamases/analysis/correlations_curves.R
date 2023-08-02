# Correlation curves of breakpoint/transposase/entropy

# Correlations to compute:
computeCorr = function(vec1, vec2){
  test = cor.test(vec1, vec2)
  return(c(as.numeric(test$estimate), test$conf.int, test$p.value))
}



data_files = Sys.glob("curves/*all-deduplicated.html.raw_data.csv")
summary_table = matrix(nrow=0, ncol=7)
i = 1
all_results = matrix(nrow=0, ncol=4)
for (file in data_files){
  gene = gsub("-curves.*", "", gsub("curves\\/", "", file))
  print(gene)
  results = read.csv(file, header=T, stringsAsFactors = F, row.names = 1)

  results.up = results[results$location=="upstream",]
  results.down = results[results$location=="downstream",]
  
  ecdf.density.up = c(gene, "upstream", "breakpoints vs.\n transposases", computeCorr(results.up$ecdf, results.up$density))
  ecdf.entropy.up = c(gene, "upstream" ,"breakpoints vs.\n block entropy", computeCorr(results.up$ecdf, results.up$entropy))
  density.entropy.up = c(gene, "upstream", "transposase vs.\n block entropy", computeCorr(results.up$density, results.up$entropy))
  ecdf.density.down = c(gene, "downstream", "breakpoints vs.\n transposases", computeCorr(results.down$ecdf, results.down$density))
  ecdf.entropy.down = c(gene, "downstream" ,"breakpoints vs.\n block entropy",computeCorr(results.down$ecdf, results.down$entropy))
  density.entropy.down = c(gene, "downstream", "transposase vs.\n block entropy",computeCorr(results.down$density, results.down$entropy))
  
  summary_table = rbind(summary_table, ecdf.density.up, ecdf.entropy.up,
        density.entropy.up, ecdf.density.down, ecdf.entropy.down, density.entropy.down)
  #results = cbind(rep(gene, nrow(results)), results )
  #all_results = rbind(all_results, results)
  #summary_table[i,] = c(gene, results)
  i = i + 1 
}

summary_table = data.frame(summary_table)
colnames(summary_table) = c("gene", "location", "corr", "estimate", "lower", "upper", "p.value")
summary_table$estimate = as.numeric(summary_table$estimate)
summary_table$lower = as.numeric(summary_table$lower)
summary_table$upper = as.numeric(summary_table$upper)
summary_table$p.value = as.numeric(summary_table$p.value)

ggplot(summary_table, aes(location, y=estimate))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_point()+
  facet_grid(corr~gene)+
  theme_basic()+
  theme(axis.text.x=element_text(angle=45, hjust=1))

p.correlations = ggplot(summary_table, aes(corr, y=estimate, group=gene, colour=gene))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_point()+
  geom_line()+
  facet_wrap(~location)+
  theme_basic()+
  ylab("correlation")+
  xlab("")
ggsave(p.correlations, filename='correlations-test-plot.pdf')
