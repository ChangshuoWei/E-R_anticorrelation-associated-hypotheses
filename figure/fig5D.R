rabies <- read.table("HRSV.txt",sep = "\t",header=T)
library(ggplot2)
library(grid)
plot(log2(rabies$expression),log2(rabies$erate))
cor.test(log2(rabies$expression),log2(rabies$erate),method = "s")
fig5d <- ggplot(rabies,aes(log2(expression),log2(erate)))
fig5d <- fig5d + 
  geom_point(color = "#de404c")+
  my_theme2
fig5d
pdf("fig5d.pdf",useDingbats = F,width = 3.6,height = 3.6)  
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout <- function(x,y)
  viewport(layout.pos.row = x,layout.pos.col = y)
print(fig5d,vp=vplayout(1,1))
dev.off()
