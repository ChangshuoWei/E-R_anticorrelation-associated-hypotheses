vaccinia2_cross <- read.table("~/E_R_anticorrelation/cross_speice/vaccinia_five/vaccinia_tree_process.txt",header = T,sep = "\t")
vaccinia2_cross_E_R <- merge(vaccinia2_cross,vaccinia_E_R,by.x = "protein",by.y = "gene")
#vaccinia2_cross_E_R <- vaccinia2_cross_E_R[vaccinia2_cross_E_R$dnds <10,]
plot(log(vaccinia2_cross_E_R$exp),log(vaccinia2_cross_E_R$dn))
cor.test(log(vaccinia2_cross_E_R$exp),log(vaccinia2_cross_E_R$dn),method = 's')
hist(vaccinia2_cross_E_R$ts.tv,breaks = 5000,xlim = c(0,))
library(ggplot2)
library(grid)
fig6c <- ggplot(vaccinia2_cross_E_R,aes(log2(exp+1),log2(dnds)))
fig6c <- fig6c + 
  geom_point(color = "#de404c")+
  my_theme2
fig6c
pdf("fig6c_2.pdf",useDingbats = F,width = 3.6,height = 3.6)  
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout <- function(x,y)
  viewport(layout.pos.row = x,layout.pos.col = y)
print(fig6c,vp=vplayout(1,1))
dev.off()
setwd("/Dell/Dell4/weics/E_R_anticorrelation/orthopoxvirus/poxvirus_project")
figS3c <- ggplot(vaccinia2_cross_E_R,aes(log2(exp+1),log2(dn)))
figS3c <- figS3c + 
  geom_point(color = "#de404c")+
  my_theme2
figS3c
pdf("figS3c.pdf",useDingbats = F,width = 3.6,height = 3.6)  
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout <- function(x,y)
  viewport(layout.pos.row = x,layout.pos.col = y)
print(fig6c,vp=vplayout(1,1))
dev.off()
cor.test(log(vaccinia2_cross_E_R$exp),log(vaccinia2_cross_E_R$dn),method = 's')
