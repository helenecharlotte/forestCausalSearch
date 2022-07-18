### rankingPerformance.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 15 2022 (13:56) 
## Version: 
## Last-Updated: Jun 23 2022 (08:10) 
##           By: Thomas Alexander Gerds
##     Update #: 157
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
rankingPerformance <- function(estimate,thecause = 1){
    ## estimate <- tar_read("ESTIMATE_ATE")
    estimate = estimate[!is.na(rank)]
    ranking = estimate[theme == "ranking",.(rank = 1:10,mean = sapply(1:10,function(r)mean(rank == r))),by = c("net","intervene","n","scale.censored","A2_T2")]
    ranking[]
}

######################################################################
### rankingPerformance.R ends here
