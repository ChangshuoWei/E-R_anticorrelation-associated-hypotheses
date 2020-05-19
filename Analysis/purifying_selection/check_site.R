non_syn_once_nucl <- read.table("abstract.txt",header = F)
names(non_syn_once_nucl) <- c("pos","refer","alter")
non_syn_once_nucl <- non_syn_once_nucl[order(non_syn_once_nucl$alter),]
non_syn_once_nucl <- non_syn_once_nucl[order(non_syn_once_nucl$refer),]
non_syn_once_nucl <- non_syn_once_nucl[order(non_syn_once_nucl$pos),]
non_syn_once_nucl$judge <- 1
for (i in 1:(nrow(non_syn_once_nucl)-1)) {
  if(non_syn_once_nucl[i,1] == non_syn_once_nucl[i+1,1] & non_syn_once_nucl[i,3] == non_syn_once_nucl[i+1,3]){
    non_syn_once_nucl[i+1,4] <- 0
  }
}
j <- 1
non_syn_once_nucl$count <- 0
while(j<nrow(non_syn_once_nucl)) {
  if(non_syn_once_nucl[j,4] ==1){
    count = 1
    h =j+1
    while(non_syn_once_nucl[h,4] == 0){
      count =count +1
      h = h+1
    }
    non_syn_once_nucl[j,5] = count
    j = j+count
  }
  else{
    j = j+1
  }
}
non_syn_once_nucl_filter <- non_syn_once_nucl[non_syn_once_nucl$count > 2,]
write.table(non_syn_once_nucl_filter,"check_var_once_test.txt",sep = "\t",quote = F,row.names = F)
system(" python main_once.py")
non_syn_once_test <- read.csv("non_syn_statistc_once_test.txt",sep = "\t")
kaks_site <- read.table('KaKs_site_paml.txt' , sep= '\t')
non_syn_once_test <- cbind(non_syn_once_test,kaks_site)
  non_syn_once_test$odd_ratio <- NA
  non_syn_once_test$p_value <- NA
for (i in 1:9) {
    non_syn_once_test[i,7] <-fisher.test(matrix(c(non_syn_once_test[i,2],non_syn_once_test[i,3],(non_syn_once_test[i,4]-non_syn_once_test[i,2]),(non_syn_once_test[i,5]-non_syn_once_test[i,3])),nrow=2))$p.value
    non_syn_once_test[i,6] <-fisher.test(matrix(c(non_syn_once_test[i,2],non_syn_once_test[i,3],(non_syn_once_test[i,4]-non_syn_once_test[i,2]),(non_syn_once_test[i,5]-non_syn_once_test[i,3])),nrow=2))$estimate
  }
write.table(non_syn_once_test,'purifying_selection_in_human.txt',sep = "\t",row.names = F,quote = F)

