data_summary<-function(data,varname,grps){
  require(plyr)
  summary_func<-function(x,col){
    c(mean=mean(x[[col]], na.rm=TRUE),
      se = sd(x[[col]], na.rm=TRUE)/length(x[[col]])^0.5)
  }
  data_sum<-ddply(data,grps,.fun=summary_func,varname)
  data_sum<-rename(data_sum,c("mean"=varname))
  return(data_sum)
}

table1 <-read.table("statistic.txt",header = T)
sars_length <- read.table("sars_length",sep= ",",header = T)
table1$length <- sars_length$length
table1$e_rate <- table1$pairwise_new/(table1$length/3 - 1)
virus_merge <- read.table("virus.merge.results",header = T)
vero_infected <- virus_merge[1:27985,]
AGM_count <- read.table("AGM_count",header = T)
names(AGM_count) <- c("transcript_id","coverage")
colnames(vero_infected)[2] <- 'transcript_id'
AGM_count <- merge(x=vero_infected,y=AGM_count)
AGM_count <- AGM_count[-10,]
AGM_evolution <- read.table("vevlet.txt",header = T,sep="\t")
plot(AGM_evolution$Macaque,AGM_evolution$query)
AGM_evolution_filter <- AGM_evolution[abs(AGM_evolution$Macaque - AGM_evolution$query) < 5,]
plot(AGM_evolution_filter$Macaque,AGM_evolution_filter$query)
cor.test(AGM_evolution_filter$Macaque,AGM_evolution_filter$query)
AGM_E_R <- merge(AGM_count,AGM_evolution_filter)
AGM_E_R$e_rate <- ((100 - AGM_E_R$Macaque)+(100 - AGM_E_R$query))/2
AGM_E_R<- AGM_E_R[AGM_E_R$e_rate < 20,]
AGM_E_R$log_coverage <- log2(AGM_E_R$coverage)
AGM_E_R$log_coverage_bin <- NA
for (i in 1:nrow(AGM_E_R)) {
  if(AGM_E_R[i,15] < 1){
    AGM_E_R[i,16] <- 1
  }
  else if(AGM_E_R[i,15] > 10){
    AGM_E_R[i,16] <- 10
  }
  else{
    AGM_E_R[i,16] <- ceiling(AGM_E_R[i,15])
  }
}
AGM_E_R$log_e_rate <- log2(AGM_E_R$e_rate+0.01)
library(grid)
library(ggplot2)
AGM_E_R_bin <- data_summary(AGM_E_R,"log_e_rate","log_coverage_bin")
plot(AGM_E_R_bin$log_e_rate)
AGM_E_R_bin$name <- "AGM"
AGM_E_R_bin$type <- "vero"
AGM_E_R_bin <- AGM_E_R_bin[,c(4,1,2,3,5)]
table1$log_leader <- log2(table1$L_reads)
table1$log_e_rate <- log2(table1$e_rate)+10
table1$log_coverage <- log(table1$nanoball_coverage)
polymorphism <- table1[,c(1,9,8)]
polymorphism$se <- 0
polymorphism$type <- "viurs"
colnames(polymorphism)[1] <- "name"
colnames(polymorphism)[2] <- "log_coverage_bin"
AGM_virus1 <- rbind(AGM_E_R_bin,polymorphism) 
g <- ggplot(AGM_virus1,aes(log_coverage_bin,log_e_rate,color = type))
fig2a <- g + 
  geom_point(size=1)+
  geom_errorbar(aes(ymin=log_e_rate-se,ymax=log_e_rate+se,width=0.2),size=0.4)+
  geom_abline(intercept = 0.37457,slope =  -0.19335,color = "#1368ac",xlim=c(0.5,15) )+
  scale_color_manual(values=c("#1368ac","#de404c"))+
  geom_text(aes(label = name))+
  my_theme2
fig2a
tpdf("fi2a.pdf",useDingbats = F,width = 4.8,height = 3.6)  
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout <- function(x,y)
  viewport(layout.pos.row = x,layout.pos.col = y)
print(fig2a,vp=vplayout(1,1))
dev.off()

cor.test(table1$log_leader,table1$log_e_rate,method = "spearman")
cor.test(AGM_E_R$log_coverage,AGM_E_R$log_e_rate,method = "spearman")
