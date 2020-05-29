fig3e <- ggplot(table1,aes(log_leader,omega))
fig3e <-fig3e + 
  geom_point(size=1,fill = "#de404c",color = "#de404c")+
  geom_text(aes(x=log_leader+0.5,label = protein,color = "#de404c"))+
  my_theme2
fig3e
pdf("fig3e.pdf",useDingbats = F,width = 4.8,height = 3.6)  
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout <- function(x,y)
  viewport(layout.pos.row = x,layout.pos.col = y)
print(fig3e,vp=vplayout(1,1))
dev.off()
cor.test(table1$log_leader,table1$omega,method = "s")
