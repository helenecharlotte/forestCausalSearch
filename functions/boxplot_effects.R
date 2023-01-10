### boxplot_effects.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 13 2022 (09:59) 
## Version: 
## Last-Updated: Jul 13 2022 (10:32) 
##           By: Thomas Alexander Gerds
##     Update #: 6
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
boxplot_effects <- function(data){
    g <- ggplot(data,aes(y = ate))
    g <- g+geom_boxplot()+theme(legend.position="none")
    g <- g+geom_hline(aes(yintercept = true.ate,color = "red"),data = data)
    g = g+ylab("Average treatment effect")+theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
    g
}


######################################################################
### boxplot_effects.R ends here
