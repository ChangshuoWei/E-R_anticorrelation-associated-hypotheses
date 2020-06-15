mumps <- read.table("mumps.txt",sep = "\t",header=T)
library(ggplot2)
library(grid)
plot(log2(mumps$expression),log2(mumps$erate))
cor.test(log2(mumps$expression),log2(mumps$erate),method = "s")
fig5f <- ggplot(mumps,aes(log2(expression),log2(erate)))
fig5f <- fig5f + 
  geom_point(color = "#de404c")+
  my_theme2
fig5f
pdf("fig5e.pdf",useDingbats = F,width = 3.6,height = 3.6)  
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout <- function(x,y)
  viewport(layout.pos.row = x,layout.pos.col = y)
print(fig5f,vp=vplayout(1,1))
dev.off()
