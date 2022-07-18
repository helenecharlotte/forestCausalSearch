### theTruth.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 10 2022 (09:18) 
## Version: 
## Last-Updated: May 20 2022 (08:45) 
##           By: Thomas Alexander Gerds
##     Update #: 89
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
plotTruth <- function(setting,
                      A1_T1,
                      A1_T2,
                      A2_T1,
                      A2_T2,
                      horizon,
                      n = 10000,
                      scale.censored){
    d = simulateData(setting = setting,
                     A1_T1 = A1_T1,
                     A1_T2 = A1_T2,
                     A2_T1 = A2_T1,
                     A2_T2 = A2_T2,
                     n = n,
                     scale.censored = scale.censored,
                     keep.latent = TRUE)
    d[,dummy := rep(1,n)]
    d[,T_treated_placebo := pmin(T1_treated_placebo,T2_treated_placebo)]
    d[,T_placebo_placebo := pmin(T1_placebo_placebo,T2_placebo_placebo)]
    p1 = with(d,ggplotify::as.grob(function(){
        plot(prodlim(Hist(T1_placebo_placebo,dummy)~1,conf.int = FALSE),xlim = c(0,horizon),type = "risk",plot.main = "Effect A1 on latent cause 1")
        legend(x = "topright",legend = c(0,1),title = "A1",col = c(1,"#E69F00"),lwd = c(2,2),cex = 1.5,bty = "n")
        plot(prodlim(Hist(T1_treated_placebo,dummy)~1,conf.int = FALSE),add = TRUE,xlim = c(0,horizon),col = "#E69F00",type = "risk")
    }))
    p2 = with(d,ggplotify::as.grob(function(){
        plot(prodlim(Hist(T_placebo_placebo,dummy)~1,conf.int = FALSE),xlim = c(0,horizon),type = "risk",plot.main = "Effect A1 on absolute risk cause 1")
        legend(x = "topright",legend = c(0,1),title = "A1",col = c(1,"#E69F00"),lwd = c(2,2),cex = 1.5,bty = "n")
        plot(prodlim(Hist(T_treated_placebo,dummy)~1,conf.int = FALSE),add = TRUE,xlim = c(0,horizon),col = "#E69F00",type = "risk")
    }))
    p3 = with(d,ggplotify::as.grob(function(){plot(prodlim(Hist(time,event)~A1,conf.int = FALSE),xlim = c(0,horizon),plot.main = "Aalen-Johansen A1")}))
    p4 = with(d,ggplotify::as.grob(function(){plot(prodlim(Hist(time,event)~A2,conf.int = FALSE),xlim = c(0,horizon),plot.main = "Aalen-Johansen A2")}))
    p5 = with(d,ggplotify::as.grob(function(){plot(prodlim(Hist(time,event)~A1,conf.int = FALSE,reverse = TRUE),xlim = c(0,horizon),plot.main = "Reverse Kaplan-Meier A1")}))
    d[,X6 := round(X6,2)]
    p6 = with(d,ggplotify::as.grob(function(){plot(prodlim(Hist(time,event)~X6,conf.int = FALSE,reverse = TRUE),xlim = c(0,horizon),plot.main = "Reverse Kaplan-Meier X6")}))
    cowplot::plot_grid(p1,p2,p3,p4,p5,p6,ncol = 2)
}

######################################################################
### theTruth.R ends here
