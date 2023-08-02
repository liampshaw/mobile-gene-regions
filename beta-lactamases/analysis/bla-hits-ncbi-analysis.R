library(dplyr)
# 22 June 2023
d = read.csv('bla-hits-ncbi.csv', header=T, stringsAsFactors = F)

d$type = gsub("-[0-9].*$", "", d$element_symbol)


summary.data = d %>% group_by(type) %>%
  summarise(n=sum(num_found)) %>%
  arrange(n, decreasing=TRUE)

print(tail(summary.data, n=50), n=50)
