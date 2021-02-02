cor.test(log2(sarsVSsars2$e_rate),log2(sarsVSsars2$pairwise_new ),method = 's')
library(ggplot2)
library(grid)
fig6a_new <- ggplot(sarsVSsars2,aes(log2(Leader.containing_reads),log2(dn.y)))
fig6a_new <- fig6a_new + 
  geom_point(color = "#de404c")+
  my_theme2
fig6a_new
pdf("figs_dn_tree.pdf",useDingbats = F,width = 3.6,height = 3.6)  
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout <- function(x,y)
  viewport(layout.pos.row = x,layout.pos.col = y)
print(fig6a_new,vp=vplayout(1,1))
dev.off()
