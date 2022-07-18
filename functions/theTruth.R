### theTruth.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 10 2022 (09:18) 
## Version: 
## Last-Updated: Jul 11 2022 (10:22) 
##           By: Thomas Alexander Gerds
##     Update #: 78
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
theTruth <- function(setting,
                     A1_T1,
                     A1_T2,
                     A2_T1,
                     A2_T2,
                     horizon,
                     n = 10000,
                     scale.censored,
                     B = 1,
                     cores = 1){
    out <- do.call("rbind",mclapply(1:B,function(b){
        d = simulateData(setting = setting,
                         A1_T1 = A1_T1,
                         A1_T2 = A1_T2,
                         A2_T1 = A2_T1,
                         A2_T2 = A2_T2,
                         n = n,
                         scale.censored = scale.censored,
                         keep.latent = TRUE)
        # average treatment effects in both worlds
        x <- d[,data.table::data.table(intervene = c("A1","A2","A1","A2","A1","A2","A1","A2"),
                                       cause = c(1,1,1,1,2,2,2,2),
                                       net = c(1,1,0,0,1,1,0,0),
                                       ate = c(
                                           # cause 1: net
                                           mean((T1_treated_A1 <= horizon) - (T1_placebo_A1 <= horizon)),
                                           mean((T1_treated_A2 <= horizon) - (T1_placebo_A2 <= horizon)),
                                           # cause 1: crude
                                           mean(T1_treated_A1 <= pmin(T2_treated_A1,horizon)) - mean(T1_placebo_A1 <= pmin(T2_placebo_A1, horizon)),
                                           mean(T1_treated_A2 <= pmin(T2_treated_A2,horizon)) - mean(T1_placebo_A2 <= pmin(T2_placebo_A2, horizon)),
                                           # cause 2: net
                                           mean((T2_treated_A1 <= horizon) - (T2_placebo_A1 <= horizon)),
                                           mean((T2_placebo_A2 <= horizon) - (T2_placebo_A2 <= horizon)),
                                           # cause 2: crude
                                           mean(T2_treated_A1 <= pmin(T1_treated_A1,horizon)) - mean(T2_placebo_A1 <= pmin(T1_placebo_A1, horizon)),
                                           mean(T2_treated_A2 <= pmin(T1_treated_A2,horizon)) - mean(T2_placebo_A2 <= pmin(T1_placebo_A2, horizon))),
                                       # calculate percent censored before tau
                                       scale.censored = scale.censored,
                                       censored.tau = round(100*mean(time <= horizon & event == 0),1))]
        ## rm(d)
        ## gc()
        x[,rep := b]
    },mc.cores = cores))
    gc()
    out = out[,data.table::data.table(ate = mean(ate),censored.tau = mean(censored.tau)),by = c("intervene","cause","net","scale.censored")]
    out[]
}

######################################################################
### theTruth.R ends here
