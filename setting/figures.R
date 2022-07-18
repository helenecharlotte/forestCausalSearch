### figures.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jun 24 2022 (14:02) 
## Version: 
## Last-Updated: Jun 24 2022 (15:04) 
##           By: Thomas Alexander Gerds
##     Update #: 4
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
try(setwd("~/research/SoftWare/grfCausalSearch/"),silent=TRUE)
library(targets)
library(tarchetypes)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggplotify)
library(cowplot)

pdf("output/estimation-performance.pdf")
b=tar_read(BOXPLOTS)
b[[5]]
dev.off()

pdf("output/ranking-performance.pdf")
ran <- tar_read(RANKING)[A2_T2%in%c(0.2,1,2)&scale.censored==0.025]
ran[,net:=factor(net,levels=c(0,1),labels=c("Crude","Net"))]
gnet=ggplot(ran[net=="Net"&intervene%in%c("A1","A2","A3")&rank==1],aes(x=n,y=mean,linetype=intervene,group=intervene))+geom_line()+geom_point()+facet_grid(~A2_T2)+ylim(c(0,1))+ylab("Frequency of rank 1")
gcrude=ggplot(ran[net=="Crude"&intervene%in%c("A1","A2","A3")&rank==1],aes(x=n,y=mean,linetype=intervene,group=intervene))+geom_line()+geom_point()+facet_grid(~A2_T2)+ylim(c(0,1))+ylab("Frequency of rank 1")
cowplot::plot_grid(gcrude+ggtitle("Crude effects"),gnet+ggtitle("Net effects"),ncol = 1)
dev.off()

######################################################################
### figures.R ends here
