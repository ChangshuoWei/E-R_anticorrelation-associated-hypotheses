#! /usr/bin/env Rscript
# -*- coding: utf-8 -*-
# @Author  : Yanming Chen
# @Version : 1.0.0


library('tidyverse')
library('data.table')

sc2.vero <- fread('out_t.protein_erate.txt', header = T, sep = '\t',
                  stringsAsFactors = F) %>% as.data.frame()
sc2.vero.plot <- sc2.vero[sc2.vero$orf != 'ORF10',]


plot.protlevel.erate <- function(indata, intag, inerate){
  colnames(indata)[grep(inerate, colnames(indata))] <- 'erate'
  indata <- na.omit(indata)
  # print(indata)
  indata$tag <- ifelse(grepl('nsp', indata$orf), 1, 0)
  outcor <- cor.test(indata$prot, indata$erate, method='s')
  # outcor <- cor.test(log(indata$prot), log(indata$erate))
  
  return(
    ggplot(indata, aes(x=log2(prot), y=log2(erate)))+
      geom_point()+
      # geom_point(aes(color=as.factor(tag)))+
      theme_classic()+
      # geom_point(data=data.frame(x=-1, y=inlow), aes(x=x, y=y), color='red', inherit.aes=F)+
      # scale_y_continuous(breaks=seq(ceiling(inlow),inup, 1))+
      xlab('log2(Protein level)')+
      ylab(sprintf('log2(Evolutionary rate in human)', inerate))+
      theme(axis.ticks=element_line(color='black'),
            axis.title=element_text(size=8, color='black'),
            axis.text=element_text(size=8, color='black'))+
      # annotate('text', label=sprintf('rho = %.2f\nP = %.4f', outcor$estimate, outcor$p.value),
      annotate('text', label=sprintf('rho = %.4f\nP = %.4f', outcor$estimate, outcor$p.value),
               x=Inf, y=-Inf, vjust=-1, hjust=1)#+
      # labs(caption=sprintf('(%s)', intag))
  )
}

{
  pdf('out3_protein_vs_erate.pdf', width=2.5, height=2.4, useDingbats=F)
  plot(plot.protlevel.erate(na.omit(sc2.vero.plot[sc2.vero.plot$orf != 'orf1ab',]), 'nsp', 'polyrate_wei'))
  dev.off()
}

## the end####




























