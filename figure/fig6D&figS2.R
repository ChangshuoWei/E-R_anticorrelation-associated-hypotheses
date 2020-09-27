healthy1_expression <- read.table("~/HERV/seq_data/healthy_1/telescope-telescope_report.tsv",sep = "\t",header = T)
healthy1_expression <- healthy1_expression[healthy1_expression$transcript != "-",]
plot(healthy_check$final_conf.x,healthy_check$final_conf.y)
cor.test(healthy_check$final_conf.x,healthy_check$final_conf.y)
healthy2_expression <- read.table("~/HERV/seq_data/healthy_2/telescope-telescope_report.tsv",sep = "\t",header = T)
healthy2_expression <- healthy2_expression[healthy2_expression$transcript != "-",]
healthy_check <- merge(healthy1_expression,healthy2_expression,by= "transcript")
plot(log(healthy_check$final_conf.x),log(healthy_check$final_conf.y))
cor.test(healthy_check$final_conf.x,healthy_check$final_conf.y)
healthy3_expression <- read.table("~/HERV/seq_data/healthy_3/telescope-telescope_report.tsv",sep = "\t",header = T)
healthy3_expression <- healthy3_expression[healthy3_expression$transcript != "-",]
healthy_check <- merge(healthy1_expression,healthy3_expression,by= "transcript")
plot(log(healthy_check$final_conf.x),log(healthy_check$final_conf.y))
cor.test(healthy_check$final_conf.x,healthy_check$final_conf.y)
healthy4_expression <- read.table("~/HERV/seq_data/healthy_4/telescope-telescope_report.tsv",sep = "\t",header = T)
healthy4_expression <- healthy4_expression[healthy4_expression$transcript != "-",]
healthy_check <- merge(healthy1_expression,healthy4_expression,by= "transcript")
plot(log(healthy_check$final_conf.x),log(healthy_check$final_conf.y))
cor.test(healthy_check$final_conf.x,healthy_check$final_conf.y)
healthy5_expression <- read.table("~/HERV/seq_data/healthy_5/telescope-telescope_report.tsv",sep = "\t",header = T)
healthy5_expression <- healthy5_expression[healthy5_expression$transcript != "-",]
healthy_check <- merge(healthy1_expression,healthy5_expression,by= "transcript")
plot(log(healthy_check$final_conf.x),log(healthy_check$final_conf.y))
cor.test(healthy_check$final_conf.x,healthy_check$final_conf.y)
healthy6_expression <- read.table("~/HERV/seq_data/healthy_6/telescope-telescope_report.tsv",sep = "\t",header = T)
healthy6_expression <- healthy6_expression[healthy6_expression$transcript != "-",]
healthy_check <- merge(healthy5_expression,healthy6_expression,by= "transcript")
plot(log(healthy_check$final_conf.x),log(healthy_check$final_conf.y))
cor.test(healthy_check$final_conf.x,healthy_check$final_conf.y)
homolog_erate3 <- read.table("homolog_erate3.txt",sep = "\t",header = T)

erv_ER6 <- merge(homolog_erate3,healthy1_expression)
erv_ER6$erate <- 2*erv_ER6$different_aa.x /(erv_ER6$query_length.x+erv_ER6$refer_length.x)
plot(log(erv_ER6$final_conf/erv_ER6$query_length.x+0.001),log(erv_ER6$erate+0.001))
cor.test((erv_ER6$final_conf/erv_ER6$query_length.x+0.01),erv_ER6$erate,method = "spearman")
erv_ER6 <- erv_ER6[order(erv_ER6$ORF_start),]
erv_ER6 <- erv_ER6[order(erv_ER6$type),]
erv_ER6 <- erv_ER6[order(erv_ER6$chrome),]
erv_ER6$dis <- 10000
for(i in 2:nrow(erv_ER6)){
  if(erv_ER6[i,18] == erv_ER6[i-1,18] & erv_ER6[i,19] == erv_ER6[i-1,19]){
    erv_ER6[i,33] = erv_ER6[i,20] - erv_ER6[i-1,21]
  }
}
plot(log(erv_ER6[erv_ER6$dis > 100,]$final_conf/erv_ER6[erv_ER6$dis > 100,]$transcript_length+0.001),log(erv_ER6[erv_ER6$dis > 100,]$erate+0.001))
cor.test((erv_ER6[erv_ER6$dis > 100,]$final_conf/erv_ER6[erv_ER6$dis > 100,]$transcript_length+0.01),erv_ER6[erv_ER6$dis > 100,]$erate,method = "spearman")
erv_ER_filer <- erv_ER6[erv_ER6$dis > 100,]
erv_ER_filer$log_expression <- log2(erv_ER_filer$final_conf/erv_ER_filer$transcript_length+0.0001)
erv_ER_filer$log_erate <- log2(erv_ER_filer$erate+0.001)
library(ggplot2)
cor.test(erv_ER_filer$log_expression,erv_ER_filer$erate,method = 's')




erv_ER_filer <- erv_ER_filer[order(erv_ER_filer$log_expression),]
erv_ER_filer$bin = NA
bin_number = 10
bin_count = floor(nrow(erv_ER_filer)/bin_number)
for (i in 1:nrow(erv_ER_filer)) {
  if(erv_ER_filer[i,34] < -10){
    erv_ER_filer[i,36] <- -10
  }
  else if(erv_ER_filer[i,34] > -5){
    erv_ER_filer[i,36] <- -5
  }
  else{
    erv_ER_filer[i,36] <- ceiling(erv_ER_filer[i,34])
  }
}
line = lm(log_erate~log_expression,erv_ER_filer)
erv_ER_bin <- data_summary(erv_ER_filer,"log_erate","bin")
fig7_1 <- ggplot(erv_ER_bin,aes(bin,log_erate))
fig7_1 <- fig7_1 + 
  geom_point(size=1,color="#de404c")+
  geom_errorbar(aes(ymin=log_erate-se,ymax=log_erate+se,width=0.2,color="#de404c"),size=0.4,)+
  geom_abline(intercept = line$coefficients[1],slope =  line$coefficients[2],color = "#de404c")+
  xlab("log2(mRNA expression level)")+
  ylab("log2(protein evolutionary rate)")+
  annotate("text", x=-7.2,y=Inf,
           label=paste("p = ",signif(cor.test(erv_ER_filer$log_expression,erv_ER_filer$erate,method = 's')$p.value,digits=2),
                       sep=" "),vjust=3.3, size=3,color="#de404c")+
  annotate("text", x=-7.2,y=Inf,
           label=paste("pho = ",signif(cor.test(erv_ER_filer$log_expression,erv_ER_filer$erate,method = 's')$estimate,digits=2), 
                       sep=" "),vjust=1.5, size=3,color="#de404c")+
  annotate("text", x=-7.2,y=Inf,
           label=paste("N =",nrow(erv_ER_filer),sep=" "), vjust=5.1, size=3,color="#de404c")+
  mytheme_violin

fig7_1


erv_ER6 <- merge(homolog_erate3,healthy2_expression)
erv_ER6$erate <- 2*erv_ER6$different_aa.x /(erv_ER6$query_length.x+erv_ER6$refer_length.x)
plot(log(erv_ER6$final_conf/erv_ER6$query_length.x+0.001),log(erv_ER6$erate+0.001))
cor.test((erv_ER6$final_conf/erv_ER6$query_length.x+0.01),erv_ER6$erate,method = "spearman")
erv_ER6 <- erv_ER6[order(erv_ER6$ORF_start),]
erv_ER6 <- erv_ER6[order(erv_ER6$type),]
erv_ER6 <- erv_ER6[order(erv_ER6$chrome),]
erv_ER6$dis <- 10000
for(i in 2:nrow(erv_ER6)){
  if(erv_ER6[i,18] == erv_ER6[i-1,18] & erv_ER6[i,19] == erv_ER6[i-1,19]){
    erv_ER6[i,33] = erv_ER6[i,20] - erv_ER6[i-1,21]
  }
}
plot(log(erv_ER6[erv_ER6$dis > 100,]$final_conf/erv_ER6[erv_ER6$dis > 100,]$transcript_length+0.001),log(erv_ER6[erv_ER6$dis > 100,]$erate+0.001))
cor.test((erv_ER6[erv_ER6$dis > 100,]$final_conf/erv_ER6[erv_ER6$dis > 100,]$transcript_length+0.01),erv_ER6[erv_ER6$dis > 100,]$erate,method = "spearman")
erv_ER_filer <- erv_ER6[erv_ER6$dis > 100,]
erv_ER_filer$log_expression <- log2(erv_ER_filer$final_conf/erv_ER_filer$transcript_length+0.0001)
erv_ER_filer$log_erate <- log2(erv_ER_filer$erate+0.001)
library(ggplot2)
cor.test(erv_ER_filer$log_expression,erv_ER_filer$erate,method = 's')




erv_ER_filer <- erv_ER_filer[order(erv_ER_filer$log_expression),]
erv_ER_filer$bin = NA
bin_number = 10
bin_count = floor(nrow(erv_ER_filer)/bin_number)
for (i in 1:nrow(erv_ER_filer)) {
  if(erv_ER_filer[i,34] < -10){
    erv_ER_filer[i,36] <- -10
  }
  else if(erv_ER_filer[i,34] > -5){
    erv_ER_filer[i,36] <- -5
  }
  else{
    erv_ER_filer[i,36] <- ceiling(erv_ER_filer[i,34])
  }
}
line = lm(log_erate~log_expression,erv_ER_filer)
erv_ER_bin <- data_summary(erv_ER_filer,"log_erate","bin")
fig7_2 <- ggplot(erv_ER_bin,aes(bin,log_erate))
fig7_2 <- fig7_2 + 
  geom_point(size=1,color="#de404c")+
  geom_errorbar(aes(ymin=log_erate-se,ymax=log_erate+se,width=0.2,color="#de404c"),size=0.4,)+
  geom_abline(intercept = line$coefficients[1],slope =  line$coefficients[2],color = "#de404c")+
  xlab("log2(mRNA expression level)")+
  ylab("log2(protein evolutionary rate)")+
  annotate("text", x=-7.2,y=Inf,
           label=paste("p = ",signif(cor.test(erv_ER_filer$log_expression,erv_ER_filer$erate,method = 's')$p.value,digits=2),
                       sep=" "),vjust=3.3, size=3,color="#de404c")+
  annotate("text", x=-7.2,y=Inf,
           label=paste("pho = ",signif(cor.test(erv_ER_filer$log_expression,erv_ER_filer$erate,method = 's')$estimate,digits=2), 
                       sep=" "),vjust=1.5, size=3,color="#de404c")+
  annotate("text", x=-7.2,y=Inf,
           label=paste("N =",nrow(erv_ER_filer),sep=" "), vjust=5.1, size=3,color="#de404c")+
  mytheme_violin

fig7_2


erv_ER6 <- merge(homolog_erate3,healthy3_expression)
erv_ER6$erate <- 2*erv_ER6$different_aa.x /(erv_ER6$query_length.x+erv_ER6$refer_length.x)
plot(log(erv_ER6$final_conf/erv_ER6$query_length.x+0.001),log(erv_ER6$erate+0.001))
cor.test((erv_ER6$final_conf/erv_ER6$query_length.x+0.01),erv_ER6$erate,method = "spearman")
erv_ER6 <- erv_ER6[order(erv_ER6$ORF_start),]
erv_ER6 <- erv_ER6[order(erv_ER6$type),]
erv_ER6 <- erv_ER6[order(erv_ER6$chrome),]
erv_ER6$dis <- 10000
for(i in 2:nrow(erv_ER6)){
  if(erv_ER6[i,18] == erv_ER6[i-1,18] & erv_ER6[i,19] == erv_ER6[i-1,19]){
    erv_ER6[i,33] = erv_ER6[i,20] - erv_ER6[i-1,21]
  }
}
plot(log(erv_ER6[erv_ER6$dis > 100,]$final_conf/erv_ER6[erv_ER6$dis > 100,]$transcript_length+0.001),log(erv_ER6[erv_ER6$dis > 100,]$erate+0.001))
cor.test((erv_ER6[erv_ER6$dis > 100,]$final_conf/erv_ER6[erv_ER6$dis > 100,]$transcript_length+0.01),erv_ER6[erv_ER6$dis > 100,]$erate,method = "spearman")
erv_ER_filer <- erv_ER6[erv_ER6$dis > 100,]
erv_ER_filer$log_expression <- log2(erv_ER_filer$final_conf/erv_ER_filer$transcript_length+0.0001)
erv_ER_filer$log_erate <- log2(erv_ER_filer$erate+0.001)
library(ggplot2)
cor.test(erv_ER_filer$log_expression,erv_ER_filer$erate,method = 's')




erv_ER_filer <- erv_ER_filer[order(erv_ER_filer$log_expression),]
erv_ER_filer$bin = NA
bin_number = 10
bin_count = floor(nrow(erv_ER_filer)/bin_number)
for (i in 1:nrow(erv_ER_filer)) {
  if(erv_ER_filer[i,34] < -10){
    erv_ER_filer[i,36] <- -10
  }
  else if(erv_ER_filer[i,34] > -5){
    erv_ER_filer[i,36] <- -5
  }
  else{
    erv_ER_filer[i,36] <- ceiling(erv_ER_filer[i,34])
  }
}
line = lm(log_erate~log_expression,erv_ER_filer)
erv_ER_bin <- data_summary(erv_ER_filer,"log_erate","bin")
fig7_3 <- ggplot(erv_ER_bin,aes(bin,log_erate))
fig7_3 <- fig7_3 + 
  geom_point(size=1,color="#de404c")+
  geom_errorbar(aes(ymin=log_erate-se,ymax=log_erate+se,width=0.2,color="#de404c"),size=0.4,)+
  geom_abline(intercept = line$coefficients[1],slope =  line$coefficients[2],color = "#de404c")+
  xlab("log2(mRNA expression level)")+
  ylab("log2(protein evolutionary rate)")+
  annotate("text", x=-7.2,y=Inf,
           label=paste("p = ",signif(cor.test(erv_ER_filer$log_expression,erv_ER_filer$erate,method = 's')$p.value,digits=2),
                       sep=" "),vjust=3.3, size=3,color="#de404c")+
  annotate("text", x=-7.2,y=Inf,
           label=paste("pho = ",signif(cor.test(erv_ER_filer$log_expression,erv_ER_filer$erate,method = 's')$estimate,digits=2), 
                       sep=" "),vjust=1.5, size=3,color="#de404c")+
  annotate("text", x=-7.2,y=Inf,
           label=paste("N =",nrow(erv_ER_filer),sep=" "), vjust=5.1, size=3,color="#de404c")+
  mytheme_violin

fig7_3


erv_ER6 <- merge(homolog_erate3,healthy4_expression)
erv_ER6$erate <- 2*erv_ER6$different_aa.x /(erv_ER6$query_length.x+erv_ER6$refer_length.x)
plot(log(erv_ER6$final_conf/erv_ER6$query_length.x+0.001),log(erv_ER6$erate+0.001))
cor.test((erv_ER6$final_conf/erv_ER6$query_length.x+0.01),erv_ER6$erate,method = "spearman")
erv_ER6 <- erv_ER6[order(erv_ER6$ORF_start),]
erv_ER6 <- erv_ER6[order(erv_ER6$type),]
erv_ER6 <- erv_ER6[order(erv_ER6$chrome),]
erv_ER6$dis <- 10000
for(i in 2:nrow(erv_ER6)){
  if(erv_ER6[i,18] == erv_ER6[i-1,18] & erv_ER6[i,19] == erv_ER6[i-1,19]){
    erv_ER6[i,33] = erv_ER6[i,20] - erv_ER6[i-1,21]
  }
}
plot(log(erv_ER6[erv_ER6$dis > 100,]$final_conf/erv_ER6[erv_ER6$dis > 100,]$transcript_length+0.001),log(erv_ER6[erv_ER6$dis > 100,]$erate+0.001))
cor.test((erv_ER6[erv_ER6$dis > 100,]$final_conf/erv_ER6[erv_ER6$dis > 100,]$transcript_length+0.01),erv_ER6[erv_ER6$dis > 100,]$erate,method = "spearman")
erv_ER_filer <- erv_ER6[erv_ER6$dis > 100,]
erv_ER_filer$log_expression <- log2(erv_ER_filer$final_conf/erv_ER_filer$transcript_length+0.0001)
erv_ER_filer$log_erate <- log2(erv_ER_filer$erate+0.001)
library(ggplot2)
cor.test(erv_ER_filer$log_expression,erv_ER_filer$erate,method = 's')




erv_ER_filer <- erv_ER_filer[order(erv_ER_filer$log_expression),]
erv_ER_filer$bin = NA
bin_number = 10
bin_count = floor(nrow(erv_ER_filer)/bin_number)
for (i in 1:nrow(erv_ER_filer)) {
  if(erv_ER_filer[i,34] < -10){
    erv_ER_filer[i,36] <- -10
  }
  else if(erv_ER_filer[i,34] > -5){
    erv_ER_filer[i,36] <- -5
  }
  else{
    erv_ER_filer[i,36] <- ceiling(erv_ER_filer[i,34])
  }
}
line = lm(log_erate~log_expression,erv_ER_filer)
erv_ER_bin <- data_summary(erv_ER_filer,"log_erate","bin")
fig7_4 <- ggplot(erv_ER_bin,aes(bin,log_erate))
fig7_4 <- fig7_4 + 
  geom_point(size=1,color="#de404c")+
  geom_errorbar(aes(ymin=log_erate-se,ymax=log_erate+se,width=0.2,color="#de404c"),size=0.4,)+
  geom_abline(intercept = line$coefficients[1],slope =  line$coefficients[2],color = "#de404c")+
  xlab("log2(mRNA expression level)")+
  ylab("log2(protein evolutionary rate)")+
  annotate("text", x=-7.2,y=Inf,
           label=paste("p = ",signif(cor.test(erv_ER_filer$log_expression,erv_ER_filer$erate,method = 's')$p.value,digits=2),
                       sep=" "),vjust=3.3, size=3,color="#de404c")+
  annotate("text", x=-7.2,y=Inf,
           label=paste("pho = ",signif(cor.test(erv_ER_filer$log_expression,erv_ER_filer$erate,method = 's')$estimate,digits=2), 
                       sep=" "),vjust=1.5, size=3,color="#de404c")+
  annotate("text", x=-7.2,y=Inf,
           label=paste("N =",nrow(erv_ER_filer),sep=" "), vjust=5.1, size=3,color="#de404c")+
  mytheme_violin

fig7_4


erv_ER6 <- merge(homolog_erate3,healthy5_expression)
erv_ER6$erate <- 2*erv_ER6$different_aa.x /(erv_ER6$query_length.x+erv_ER6$refer_length.x)
plot(log(erv_ER6$final_conf/erv_ER6$query_length.x+0.001),log(erv_ER6$erate+0.001))
cor.test((erv_ER6$final_conf/erv_ER6$query_length.x+0.01),erv_ER6$erate,method = "spearman")
erv_ER6 <- erv_ER6[order(erv_ER6$ORF_start),]
erv_ER6 <- erv_ER6[order(erv_ER6$type),]
erv_ER6 <- erv_ER6[order(erv_ER6$chrome),]
erv_ER6$dis <- 10000
for(i in 2:nrow(erv_ER6)){
  if(erv_ER6[i,18] == erv_ER6[i-1,18] & erv_ER6[i,19] == erv_ER6[i-1,19]){
    erv_ER6[i,33] = erv_ER6[i,20] - erv_ER6[i-1,21]
  }
}
plot(log(erv_ER6[erv_ER6$dis > 100,]$final_conf/erv_ER6[erv_ER6$dis > 100,]$transcript_length+0.001),log(erv_ER6[erv_ER6$dis > 100,]$erate+0.001))
cor.test((erv_ER6[erv_ER6$dis > 100,]$final_conf/erv_ER6[erv_ER6$dis > 100,]$transcript_length+0.01),erv_ER6[erv_ER6$dis > 100,]$erate,method = "spearman")
erv_ER_filer <- erv_ER6[erv_ER6$dis > 100,]
erv_ER_filer$log_expression <- log2(erv_ER_filer$final_conf/erv_ER_filer$transcript_length+0.0001)
erv_ER_filer$log_erate <- log2(erv_ER_filer$erate+0.001)
library(ggplot2)
cor.test(erv_ER_filer$log_expression,erv_ER_filer$erate,method = 's')




erv_ER_filer <- erv_ER_filer[order(erv_ER_filer$log_expression),]
erv_ER_filer$bin = NA
bin_number = 10
bin_count = floor(nrow(erv_ER_filer)/bin_number)
for (i in 1:nrow(erv_ER_filer)) {
  if(erv_ER_filer[i,34] < -10){
    erv_ER_filer[i,36] <- -10
  }
  else if(erv_ER_filer[i,34] > -5){
    erv_ER_filer[i,36] <- -5
  }
  else{
    erv_ER_filer[i,36] <- ceiling(erv_ER_filer[i,34])
  }
}
line = lm(log_erate~log_expression,erv_ER_filer)
erv_ER_bin <- data_summary(erv_ER_filer,"log_erate","bin")
fig7_5 <- ggplot(erv_ER_bin,aes(bin,log_erate))
fig7_5 <- fig7_5 + 
  geom_point(size=1,color="#de404c")+
  geom_errorbar(aes(ymin=log_erate-se,ymax=log_erate+se,width=0.2,color="#de404c"),size=0.4,)+
  geom_abline(intercept = line$coefficients[1],slope =  line$coefficients[2],color = "#de404c")+
  xlab("log2(mRNA expression level)")+
  ylab("log2(protein evolutionary rate)")+
  annotate("text", x=-7.2,y=Inf,
           label=paste("p = ",signif(cor.test(erv_ER_filer$log_expression,erv_ER_filer$erate,method = 's')$p.value,digits=2),
                       sep=" "),vjust=3.3, size=3,color="#de404c")+
  annotate("text", x=-7.2,y=Inf,
           label=paste("pho = ",signif(cor.test(erv_ER_filer$log_expression,erv_ER_filer$erate,method = 's')$estimate,digits=2), 
                       sep=" "),vjust=1.5, size=3,color="#de404c")+
  annotate("text", x=-7.2,y=Inf,
           label=paste("N =",nrow(erv_ER_filer),sep=" "), vjust=5.1, size=3,color="#de404c")+
  mytheme_violin

fig7_5


erv_ER6 <- merge(homolog_erate3,healthy6_expression)
erv_ER6$erate <- 2*erv_ER6$different_aa.x /(erv_ER6$query_length.x+erv_ER6$refer_length.x)
plot(log(erv_ER6$final_conf/erv_ER6$query_length.x+0.001),log(erv_ER6$erate+0.001))
cor.test((erv_ER6$final_conf/erv_ER6$query_length.x+0.01),erv_ER6$erate,method = "spearman")
erv_ER6 <- erv_ER6[order(erv_ER6$ORF_start),]
erv_ER6 <- erv_ER6[order(erv_ER6$type),]
erv_ER6 <- erv_ER6[order(erv_ER6$chrome),]
erv_ER6$dis <- 10000
for(i in 2:nrow(erv_ER6)){
  if(erv_ER6[i,18] == erv_ER6[i-1,18] & erv_ER6[i,19] == erv_ER6[i-1,19]){
    erv_ER6[i,33] = erv_ER6[i,20] - erv_ER6[i-1,21]
  }
}
plot(log(erv_ER6[erv_ER6$dis > 100,]$final_conf/erv_ER6[erv_ER6$dis > 100,]$transcript_length+0.001),log(erv_ER6[erv_ER6$dis > 100,]$erate+0.001))
cor.test((erv_ER6[erv_ER6$dis > 100,]$final_conf/erv_ER6[erv_ER6$dis > 100,]$transcript_length+0.01),erv_ER6[erv_ER6$dis > 100,]$erate,method = "spearman")
erv_ER_filer <- erv_ER6[erv_ER6$dis > 100,]
erv_ER_filer$log_expression <- log2(erv_ER_filer$final_conf/erv_ER_filer$transcript_length+0.0001)
erv_ER_filer$log_erate <- log2(erv_ER_filer$erate+0.001)
library(ggplot2)
cor.test(erv_ER_filer$log_expression,erv_ER_filer$erate,method = 's')




erv_ER_filer <- erv_ER_filer[order(erv_ER_filer$log_expression),]
erv_ER_filer$bin = NA
bin_number = 10
bin_count = floor(nrow(erv_ER_filer)/bin_number)
for (i in 1:nrow(erv_ER_filer)) {
  if(erv_ER_filer[i,34] < -10){
    erv_ER_filer[i,36] <- -10
  }
  else if(erv_ER_filer[i,34] > -5){
    erv_ER_filer[i,36] <- -5
  }
  else{
    erv_ER_filer[i,36] <- ceiling(erv_ER_filer[i,34])
  }
}
line = lm(log_erate~log_expression,erv_ER_filer)
erv_ER_bin <- data_summary(erv_ER_filer,"log_erate","bin")
fig7_6 <- ggplot(erv_ER_bin,aes(bin,log_erate))
fig7_6 <- fig7_6 + 
  geom_point(size=1,color="#de404c")+
  geom_errorbar(aes(ymin=log_erate-se,ymax=log_erate+se,width=0.2,color="#de404c"),size=0.4,)+
  geom_abline(intercept = line$coefficients[1],slope =  line$coefficients[2],color = "#de404c")+
  xlab("log2(mRNA expression level)")+
  ylab("log2(protein evolutionary rate)")+
  annotate("text", x=-7.2,y=Inf,
           label=paste("p = ",signif(cor.test(erv_ER_filer$log_expression,erv_ER_filer$erate,method = 's')$p.value,digits=2),
                       sep=" "),vjust=3.3, size=3,color="#de404c")+
  annotate("text", x=-7.2,y=Inf,
           label=paste("pho = ",signif(cor.test(erv_ER_filer$log_expression,erv_ER_filer$erate,method = 's')$estimate,digits=2), 
                       sep=" "),vjust=1.5, size=3,color="#de404c")+
  annotate("text", x=-7.2,y=Inf,
           label=paste("N =",nrow(erv_ER_filer),sep=" "), vjust=5.1, size=3,color="#de404c")+
  mytheme_violin

fig7_6

library("grid")
pdf("fig7_bin6.pdf",width=20,height=20,useDingbats = F)
print(fig7_1,vp=viewport(.09,.08,x=.2,y=.8))
print(fig7_2,vp=viewport(.09,.08,x=.3,y=.8))
print(fig7_3,vp=viewport(.09,.08,x=.4,y=.8))
print(fig7_4,vp=viewport(.09,.08,x=.2,y=.7))
print(fig7_5,vp=viewport(.09,.08,x=.3,y=.7))
print(fig7_6,vp=viewport(.09,.08,x=.4,y=.7))

dev.off()

