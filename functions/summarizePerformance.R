### summarizePerformance.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 15 2022 (13:56) 
## Version: 
## Last-Updated: Jul 18 2022 (08:07) 
##           By: Thomas Alexander Gerds
##     Update #: 154
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
summarizePerformance <- function(truth,estimate,variable = "A1",thecause = 1){
    ## estimate <- tar_read("ESTIMATE_ATE")
    ## truth <- tar_read("TRUTH")
    estimate = estimate[intervene == variable]
    truth = truth[intervene == variable]
    ## setnames(estimate,"time","horizon")
    setkeyv(estimate,c("net","formula","horizon","scale.censored","A1_T1","A1_T2","A2_T1","A2_T2"))
    setkeyv(truth,c("net","formula","horizon","scale.censored","A1_T1","A1_T2","A2_T1","A2_T2"))
    # merge true ate: first calculate average true.ate because,
    # e.g., there are repeated values for different scale.censored but scale.censored does not affect true ATE
    t.ate = truth[cause == thecause,.(true.ate = mean(ate)),keyby = c("net","formula","horizon","A1_T1","A1_T2","A2_T1","A2_T2")]
    # merge percentage censored before horizon
    t.cens = truth[cause == thecause,.(censored.tau = mean(censored.tau)),keyby = c("formula","horizon","scale.censored")]
    estimate = merge(estimate,t.ate,by = c("net","formula","horizon","A1_T1","A1_T2","A2_T1","A2_T2"),all.y = FALSE)
    estimate = merge(estimate,t.cens,by = c("formula","horizon","scale.censored"),all.y = FALSE)
    perf = estimate[,.(repetitions = .N,
                       true.ate = true.ate[1],
                       mean.ate = mean(ate,na.rm = TRUE),
                       bias = mean(ate-true.ate,na.rm = TRUE),
                       abs.bias = mean(abs(ate-true.ate),na.rm = TRUE),
                       root.mse = sqrt(mean(ate-true.ate,na.rm = TRUE)^2),
                       sd = sd(ate,na.rm = TRUE),
                       mean.se = mean(se,na.rm = TRUE),
                       coverage = mean(true.ate >= lower&true.ate <= upper,na.rm = TRUE),
                       NAs = sum(is.na(ate))),
                    by = c("theme","formula","horizon","net","n","method","scale.censored","num.trees","A1_T1","A1_T2","A2_T1","A2_T2")]
    x = truth[,.(censored.tau,scale.censored)]
    x = x[!duplicated(scale.censored)]
    x[,scale.censored := factor(as.character(scale.censored))]
    perf[,scale.censored := factor(as.character(scale.censored))]
    setkey(perf,scale.censored)
    setkey(x,scale.censored)
    perf = x[perf]
    perf[]
}

######################################################################
### summarizePerformance.R ends here
