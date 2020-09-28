library(ggplot2)
library(grid)
AGM_E_R_top5 <- AGM_E_R[order(AGM_E_R$coverage,decreasing = T),]
AGM_E_R_range174 <- AGM_E_R_top5[AGM_E_R_top5$coverage > 147.6,]

AGM_E_R_top5 <- AGM_E_R_top5[c(1:1915),]
for(i in 1:nrow(AGM_E_R_top5)){
  if(AGM_E_R_top5[i,15] > 7.4 & AGM_E_R_top5[i,15] <= 7.8){
    AGM_E_R_top5[i,16] <- 7.6
  }
  else if(AGM_E_R_top5[i,15] > 7.8 & AGM_E_R_top5[i,15] <= 8.2){
    AGM_E_R_top5[i,16] <- 8.0
  }
  else if(AGM_E_R_top5[i,15] > 8.2 & AGM_E_R_top5[i,15] <= 8.6){
    AGM_E_R_top5[i,16] <- 8.4
  }
  else if(AGM_E_R_top5[i,15] > 8.6 & AGM_E_R_top5[i,15] <= 9.0){
    AGM_E_R_top5[i,16] <- 8.8
  }
  else if(AGM_E_R_top5[i,15] > 9.0 & AGM_E_R_top5[i,15] <= 9.4){
    AGM_E_R_top5[i,16] <- 9.2
  }
  else if(AGM_E_R_top5[i,15] > 9.4 & AGM_E_R_top5[i,15] <= 9.8){
    AGM_E_R_top5[i,16] <- 9.6
  }
  else if(AGM_E_R_top5[i,15] > 9.8){
    AGM_E_R_top5[i,16] <- 9.6
  }
}
library(ggplot2)




AGM_E_R_top5_bin <- data_summary(AGM_E_R_top5,"log_e_rate","log_coverage_bin")
plot(AGM_E_R_top5_bin$log_coverage_bin,AGM_E_R_top5_bin$log_e_rate)
plot(AGM_E_R_top5$log_coverage,AGM_E_R_top5$log_e_rate)
cor.test(AGM_E_R_top5$coverage,AGM_E_R_top5$e_rate,method = 's')
figs1a <- ggplot(AGM_E_R_top5_bin,aes(log_coverage_bin,log_e_rate))
figs1a <- figs1a + 
  geom_point(size=1,color="#de404c") +
  geom_errorbar(aes(ymin=log_e_rate-se,ymax=log_e_rate+se,width=0.2,color="#de404c"),size=0.4,)+
  xlab("log2(mRNA expression level)")+
  ylab("log2(protein evolutionary rate)")+
  annotate("text", x=9.5,y=Inf,
           label=paste("p = ",signif(cor.test(AGM_E_R_top5$log_coverage,AGM_E_R_top5$log_e_rate,method = 's')$p.value,digits=2),
                       sep=" "),vjust=3.3, size=3,color="#de404c")+
  annotate("text", x=9.5,y=Inf,
           label=paste("pho = ",signif(cor.test(AGM_E_R_top5$log_coverage,AGM_E_R_top5$log_e_rate,method = 's')$estimate,digits=2), 
                       sep=" "),vjust=1.5, size=3,color="#de404c")+
  annotate("text", x=9.5,y=Inf,
           label=paste("N =",nrow(AGM_E_R_top5),sep=" "), vjust=5.1, size=3,color="#de404c")+
  mytheme_violin
figs1a
