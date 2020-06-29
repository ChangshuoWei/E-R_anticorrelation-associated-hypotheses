measles <- read.table("measles.txt",sep = "\t",header=T)
library(ggplot2)
library(grid)
plot(log2(measles$mRNA),log2(measles$erate))
cor.test(log2(measles$mRNA),log2(measles$erate),method = "s")
fig5g <- ggplot(measles,aes(log2(mRNA),log2(erate)))
fig5g <- fig5g + 
  geom_point(color = "#de404c")+
  my_theme2
fig5g
pdf("fig5g.pdf",useDingbats = F,width = 3.6,height = 3.6)  
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout <- function(x,y)
  viewport(layout.pos.row = x,layout.pos.col = y)
print(fig5g,vp=vplayout(1,1))
dev.off()
