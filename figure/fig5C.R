newcastle <- read.table("HPIV.txt",sep = "\t",header=T)
library(ggplot2)
library(grid)
plot(log2(newcastle$expression),log2(newcastle$erate))
cor.test(log2(newcastle$expression),log2(newcastle$erate),method = "s")
fig5e <- ggplot(newcastle,aes(log2(expression),log2(erate)))
fig5e <- fig5e + 
  geom_point(color = "#de404c")+
  my_theme2
fig5e
pdf("fig5e.pdf",useDingbats = F,width = 3.6,height = 3.6)  
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout <- function(x,y)
  viewport(layout.pos.row = x,layout.pos.col = y)
print(fig5e,vp=vplayout(1,1))
dev.off()
