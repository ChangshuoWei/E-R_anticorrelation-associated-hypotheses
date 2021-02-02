herper_expression <- read.table("/Dell/Dell4/weics/E_R_anticorrelation/herperviral/beta/GSE99823_FPKM13.txt",header = T)
herper_evolutionary <- read.table("/Dell/Dell4/weics/E_R_anticorrelation/herperviral/beta/erate_statistic_protein.txt")
names(herper_evolutionary) <- c("Gene",'poly','erate')
herper_E_R <- merge(herper_evolutionary,herper_expression,by.x = "Gene",by.y = "Gene")
cor.test(log(herper_E_R$SS.WT_2dpi),log(herper_E_R$erate),method = 's')
plot(log(herper_E_R$SS.WT_2dpi),log(herper_E_R$erate))
hist(log(herper_E_R$erate),breaks = 150)
plot(log(herper_E_R$SS.WT_2dpi),log(herper_E_R$SS.WT_6dpi))

herper_evolutionary2 <- read.table("/Dell/Dell4/weics/E_R_anticorrelation/herperviral/beta/needle_erate/erate_stat.txt")
names(herper_evolutionary2) <- c("gene_name",'erate')
name_temp <- read.table("/Dell/Dell4/weics/E_R_anticorrelation/herperviral/beta/needle_erate/name.txt",sep = "\t")
names(name_temp) <- c("gene_name",'description')
herper_evolutionary2<- merge(name_temp,herper_evolutionary2)
name_temp2 <- read.table("/Dell/Dell4/weics/E_R_anticorrelation/herperviral/beta/NC_expression_name.txt",sep = "\t")
names(name_temp2) <- c("Gene","description","pro_id","strand")
herper_evolutionary2<- merge(name_temp2,herper_evolutionary2)
herper_evolutionary_check <- merge(herper_evolutionary[,c(1,3)],herper_evolutionary2[,c(2,6)],by.x = "Gene",by.y = "Gene")
plot(log2(herper_evolutionary_check$erate.x),log2(herper_evolutionary_check$erate.y))
cor.test(herper_evolutionary_check$erate.x,herper_evolutionary_check$erate.y)


herper_E_R <- merge(herper_evolutionary2,herper_expression,by.x = "Gene",by.y = "Gene")
cor.test(log(herper_E_R$SS.WT_2dpi),log(herper_E_R$erate),method = 's')
plot(log(herper_E_R$SS.WT_2dpi),log(herper_E_R$erate))
hist(log(herper_E_R$erate),breaks = 150)
plot(log(herper_E_R$SS.WT_2dpi),log(herper_E_R$NS.WT_2dpi))
cor.test((herper_E_R$NS.WT_2dpi),(herper_E_R$SS.WT_2dpi))

fig5I <- ggplot(herper_E_R,aes(log2(SS.WT_2dpi),log2(erate)))
fig5I <- fig5I+
  geom_point(color = "#de404c")+
  my_theme2
fig5I
pdf("fig5i.pdf",useDingbats = F,width = 3.6,height = 3.6)  
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout <- function(x,y)
  viewport(layout.pos.row = x,layout.pos.col = y)
print(fig5I,vp=vplayout(1,1))
dev.off()


herper_cross <- read.table("~/E_R_anticorrelation/cross_speice/herpes/PanineVSbeta.txt",header = T,sep = "\t")
herper_cross_E_R <- merge(herper_cross,herper_E_R,by.x = "protein",by.y = "Gene")
cor.test(log(herper_cross_E_R$SS.WT_2dpi),log(herper_cross_E_R$dnds),method = 's')
plot(log(herper_cross_E_R$SS.WT_2dpi),log(herper_cross_E_R$dnds))
hist(herper_cross_E_R$ts.tv,breaks = 50)

library(ggplot2)
library(grid)

fig6b <- ggplot(herper_cross_E_R,aes(log2(SS.WT_2dpi),log2(dnds)))
fig6b <- fig6b + 
  geom_point(color = "#de404c")+
  my_theme2
fig6b
pdf("fig6b.pdf",useDingbats = F,width = 3.6,height = 3.6)  
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout <- function(x,y)
  viewport(layout.pos.row = x,layout.pos.col = y)
print(fig6b,vp=vplayout(1,1))
dev.off()
