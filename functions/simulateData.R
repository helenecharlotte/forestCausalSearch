### simulateData.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 15 2022 (20:29) 
## Version: 
## Last-Updated: Apr 29 2022 (07:19) 
##           By: Thomas Alexander Gerds
##     Update #: 28
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
simulateData <- function(setting,
                         A1_T1,
                         A1_T2,
                         A2_T1,
                         A2_T2,
                         n,
                         scale.censored,
                         keep.latent = TRUE){
    # adapt treatment effects
    setting$treatment.effect <- list("A1_T1" = A1_T1,
                                     "A1_T2" = A1_T2,
                                     "A2_T1" = A2_T1,
                                     "A2_T2" = A2_T2)
    # latent variable model
    model = do.call("lavaModel",setting)
    if (is.infinite(scale.censored)) # FIXME: achieve uncensored by setting all C values higher than 10000 (hoping that the horizon is below 10000) 
        lava::distribution(model,"C") <- lava::uniform.lvm(a = 10000,b = 10000+1)
    else
        lava::distribution(model,"C") <- lava::coxWeibull.lvm(scale = scale.censored)
    out = data.table(lava::sim(model,n = n))
    # remove latent variables
    if (keep.latent[1] == FALSE){
        out[,c("T1_placebo_placebo","T1_treated_placebo","T1_placebo_treated","T1_treated_treated","T2_placebo_placebo","T2_treated_placebo","T2_placebo_treated","T2_treated_treated","C","T1","T2") := NULL]
    }
    out[]
}

######################################################################
### simulateData.R ends here
