library(ggplot2)
library(grid)
test = AGM_E_R[floor(runif(10,1,9576)),]
cor.test(test$coverage,test$e_rate,method = "spearman")
distri <- data.frame(1:10000)
for(j in 1: 10000){
  test = AGM_E_R[floor(runif(9,1,9576)),]
  rho = cor.test(test$coverage,test$e_rate,method = "spearman")$estimate
  distri[j,1] <- rho
}
quantile(distri$X1.1000,probs = seq(0,1,0.05))
plot(density(distri$X1.1000))
sum(distri$X1.1000 >=0.6)
sum(AGM_E_R$coverage)
colnames(distri)[1] <- "rho"
fig3f <- ggplot(distri,aes(rho))
fig3f <- fig3f + 
  geom_histogram(bins = 25,color = "#1368ac",fill = '#1368ac')+
  geom_vline(xintercept = 0.6,color="#de404c")+
  my_theme2
fig3f
pdf("fig3f.pdf",useDingbats = F,width = 4.5,height = 3.2)  
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout <- function(x,y)
  viewport(layout.pos.row = x,layout.pos.col = y)
print(fig3f,vp=vplayout(1,1))
dev.off()
