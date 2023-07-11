library(ggplot2)
library(cowplot)

# Upstream
d.upstream = data.frame(dist=c(0, 0.16, 0.16, 1.34, 1.34, 1.8, 1.8, 3.38), 
               ecdf=c(0, 0, 0.6,0.6, 0.7, 0.7,1, 1))
p.up = ggplot(d.upstream, aes(dist, ecdf))+
  scale_x_reverse()+
  scale_y_reverse()+
  geom_path()+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.line = element_blank(),
        panel.border = element_blank(),
        axis.title=element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())+
  theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA))

# save as pdf with right dimensions for schematic

d.down = data.frame(dist=c(0, 0.43,0.43, 1.96, 1.96, 3.45, 3.45),
                    ecdf=c(0,0,0.6,0.6,0.7,0.7,1))

# Downstream
p.down = ggplot(d.down, aes(dist, ecdf))+scale_y_reverse()+geom_path()+theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.line = element_blank(),
        panel.border = element_blank(),
        axis.title=element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())+
  theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA))

p.combined = cowplot::plot_grid(p.up, p.down)

ggsave(p.combined, file='schematic-breakpoint.pdf', width=13, height=3,
      )


-(2/5 * log(2/5) + 3/5*log(3/5))
-(1/5 * log(1/5) + 1/5 * log(1/5) + 3/5*log(3/5))
-(1/5 * log(1/5) *5)

entropy.upstream = data.frame(dist=c(0, 0.16, 0.16, 1.34, 1.34, 1.8,1.8,3.38), 
                        entropy=c(0, 0, 0.6730117,0.6730117, 0.9502705,0.9502705, 1.609438, 1.609438))

-(1/5 * log(1/5) + 1/5 * log(1/5) + 3/5*log(3/5))

entropy.downstream = data.frame(dist=c(0, 0.43,0.43, 1.96, 1.96, 3.45),
                                entropy=c(0,0, 0.6730117, 0.6730117, 0.6730117, 0.6730117))
p.entropy.down = ggplot(entropy.downstream, aes(dist, entropy))+geom_path()+theme_bw()+
  ylim(c(0,1.7))+
  theme(panel.grid = element_blank())+
  theme(axis.line = element_blank(),
        panel.border = element_blank(),
        axis.title=element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())+
  theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA))
p.entropy.up = ggplot(entropy.upstream, aes(dist, entropy))+geom_path()+theme_bw()+
  scale_x_reverse()+
  ylim(c(0,1.7))+
  theme(panel.grid = element_blank())+
  theme(axis.line = element_blank(),
        panel.border = element_blank(),
        axis.title=element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())+
  theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA))
p.combined.entropy = cowplot::plot_grid(p.entropy.up, p.entropy.down)

ggsave(p.combined.entropy, file='schematic-entropy.pdf', width=13, height=3,
)
