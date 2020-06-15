mers <- read.table("mers.txt",sep = "\t",header=T)
cor.test(mers$junction3,mers$e_rate,method = "s")
library(ggplot2)
library(grid)


fig5a <- ggplot(mers,aes(log2(junction3),log2(e_rate)))
fig5a <- fig5a + 
  geom_point(color = "#de404c")+
  my_theme2
fig5a
pdf("fig5a.pdf",useDingbats = F,width = 3.6,height = 3.6)  
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout <- function(x,y)
  viewport(layout.pos.row = x,layout.pos.col = y)
print(fig5a,vp=vplayout(1,1))
dev.off()
