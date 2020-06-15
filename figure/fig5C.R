ebola <- read.table("ebola.txt",header = T)
library(ggplot2)
library(grid)
fig5c <- ggplot(ebola,aes(log2(TPM),log2(erate)))
fig5c <- fig5c + 
  geom_point(color = "#de404c")+
  my_theme2
fig5c
pdf("fig5c.pdf",useDingbats = F,width = 3.6,height = 3.6)  
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout <- function(x,y)
  viewport(layout.pos.row = x,layout.pos.col = y)
print(fig5c,vp=vplayout(1,1))
dev.off()
cor.test(ebola$TPM,ebola$erate,method = 's')
