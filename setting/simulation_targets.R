### simulation_targets.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: May  5 2022 (11:08) 
## Version: 
## Last-Updated: Jan 10 2023 (14:31) 
##           By: Thomas Alexander Gerds
##     Update #: 201
#----------------------------------------------------------------------
## 
### Commentary:
# 
# effects on crude and net probabilities 
# independent censoring
# varying sample size
# varying percentage of censored (before horizon)
# varying treatment effects
# ranking performance
# 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(data.table)
## MCCORES <- 1
## MCCORES <- 5
MCCORES <- 25
setDTthreads(1)
fixed <- list(event.times = c("T1","T2"),
              treatments = c("A1" = .4,"A2" = .3, "A3" = .3,"A4" = .4,"A5" = .5,"A6" = .2,"A7" = .7,"A8" = .8,"A9" = .9,"A10" = .1),
              binary.covariates = c("X1" = .1,"X2" = .2,"X3" = .3,"X4" = .4,"X5" = .5),
              normal.covariates = paste0("X",6:7),
              quadratic.covariates = "X6")
fixed_target <- tar_target(FIXED,fixed,deployment = "main")
# basic setting: no covariate effects on censoring time 
formula1 <- list(T1 ~ f(A4,.3) + f(A5,.7) + f(X1,1) + f(X2,.3) + f(X6,-.5),
                 T2  ~ f(A5, -.3)+ f(X1,-.1) + f(X2,.6) + f(X6,.1),
                 C ~ f(A1,0),
                 A1 ~ f(X1,-1) + f(X6,.7) + f(A7,.2))
# no other treatment effects on T1 or T2 
formula_ranking <- list(T1 ~ f(X1,1) + f(X2,.3) + f(X6,.3),
                        T2  ~ f(X1,-.1) + f(X2,.6) + f(X6,.1),
                        C ~ f(A1,0),
                        A1 ~ f(X1,-1) + f(X6,.7) + f(A7,.2))
# covariate dependent censoring
formula_cens <- list(T1 ~ f(A4,.3) + f(A5,.7) + f(X1,1) + f(X2,.3) + f(X6,-.5),
                     T2  ~ f(A5, -.3)+ f(X1,-.1) + f(X2,.6) + f(X6,.1),
                     C ~ f(A1,.1)+f(A3,-0.3) + f(X2,-0.4)+f(X6,.1),
                     A1 ~ f(X1,-1) + f(X6,.7) + f(A7,.2))
# misspecified parametric models
formula_misspecified =  list(T1 ~ f(A4,.3) + f(A5,.7) + f(X1,1) + f(X2,.3) + f(X6,.3)+f(X6_2,.4),
                             T2  ~ f(A5, -.3)+ f(X1,-.1) + f(X2,.6) + f(X6,.1),
                             C ~ f(A1,0),
                             A1 ~ f(X1,-1) + f(X6,.7) + f(A7,.2))
formula_settings <- list("formula1" = formula1,
                         "formula_ranking" = formula_ranking,
                         "formula_misspecified" = formula_misspecified,
                         "formula_cens" = formula_cens)
# crude treatment effect
varying_crude <- data.table::CJ(A1_T1 = 1.25,
                                A1_T2 = c(.8,1,1.25),
                                A2_T1 = 1,
                                A2_T2 = 1.25,
                                scale.censored = 1/40,
                                sample.size = c(500,1000,2000,5000),
                                horizon = 5,
                                setting = "formula1",
                                method = "causal_forest",
                                weighter = "ranger",
                                net = FALSE,
                                treat = "A1",
                                num.trees = 50)
# net effects 
varying_net <- data.table::CJ(A1_T1 = c(.8,1,1.25),
                              A1_T2 = 0.8,
                              A2_T1 = 1,
                              A2_T2 = 1.25,
                              scale.censored = 1/40,
                              sample.size = 5000,
                              horizon = 5,
                              setting = "formula1",
                              method = "causal_forest",
                              weighter = "ranger",
                              net = c(FALSE,TRUE),
                              treat = "all",
                              num.trees = 50)
# effect censored
varying_censored <- data.table::CJ(A1_T1 = 1.25,
                                   A1_T2 = 1,
                                   A2_T1 = 1,
                                   A2_T2 = 1.25,
                                   scale.censored = c(-Inf,1/40,1/25),
                                   sample.size = c(500,1000,2000,5000),
                                   horizon = 2.5,
                                   setting = c("formula1","formula_cens"),
                                   method = "causal_forest",
                                   weighter = "ranger",
                                   net = FALSE,
                                   treat = "A1",
                                   num.trees = 50)
# misspecified parametric models
varying_misspecified <- data.table::CJ(A1_T1 = 1.25,
                                       A1_T2 = 1,
                                       A2_T1 = 1,
                                       A2_T2 = 1.25,
                                       scale.censored = 1/40,
                                       ## sample.size = 5000,
                                       sample.size = c(500,1000,2000,5000),
                                       horizon = 5,
                                       setting = "formula_misspecified",
                                       method = c("causal_forest","CSC","FGR"),
                                       weighter = "ranger",
                                       net = FALSE,
                                       treat = "A1",
                                       num.trees = 50)
varying_notmisspecified <- data.table::CJ(A1_T1 = 1.25,
                                       A1_T2 = 1,
                                       A2_T1 = 1,
                                       A2_T2 = 1.25,
                                       scale.censored = 1/40,
                                       sample.size = c(500,1000,2000,5000),
                                       horizon = 5,
                                       setting = "formula1",
                                       method = c("causal_forest","CSC","FGR"),
                                       weighter = "ranger",
                                       net = FALSE,
                                       treat = "A1",
                                       num.trees = 50)
# ranking treatments
varying_ranking <- data.table::CJ(A1_T1 = 1.25,
                                  A1_T2 = 1,
                                  A2_T1 = 1,
                                  A2_T2 = c(.2,.8,1,1.25,2),
                                  scale.censored = c(-Inf,1/40),
                                  sample.size = c(500,1000,2000,5000),
                                  horizon = 5,
                                  setting = "formula_ranking",
                                  method = "causal_forest",
                                  weighter = "ranger",
                                  net = c(FALSE,TRUE),
                                  treat = "all",
                                  num.trees = 100)
varying <-  rbindlist(list(varying_crude[,theme := "crude_effect"],
                           varying_net[,theme := "net_effect"],
                           varying_censored[,theme := "censoring"],
                           varying_misspecified[,theme := "misspecified"],
                           varying_notmisspecified[,theme := "not_misspecified"],
                           varying_ranking[,theme := "ranking"]))
varying[sample.size == 5000,num.trees := 50]
varying_target <- tar_target(VARYING,
                             varying,
                             deployment = "main")
varying[,reps := 2002]
varying[sample.size ==  5000,reps := 1001]

# ---------------------------------------------------------------------
# Calculation of true ATE values
# ---------------------------------------------------------------------
#   crude: real world with competing risks
#     net: hypothetical world where competing risks are eliminated
truth_varying <- tar_map(
    # truth is affected by % censored and true effect
    values = unique(varying[,.(scale.censored,A1_T1,A1_T2,A2_T1,A2_T2,setting,horizon,theme)]),
    tar_target(VARYING_TRUTH,{
        ## lavaModel;
        ## simulateData;
        set <- FIXED
        set$formula.list = formula_settings[[setting]]
        y = theTruth(setting = set,
                     A1_T1 = A1_T1,
                     A1_T2 = A1_T2,
                     A2_T1 = A2_T1,
                     A2_T2 = A2_T2,
                     scale.censored = scale.censored,
                     horizon = horizon,
                     B = 50,
                     cores = MCCORES,
                     n = 100000)
        y = cbind(y,
                  A1_T1 = A1_T1,
                  A1_T2 = A1_T2,
                  A2_T1 = A2_T1,
                  A2_T2 = A2_T2,
                  horizon = horizon,
                  theme = theme,
                  formula = setting)
        y},
        deployment = "main"))

truth <- tar_combine(TRUTH, truth_varying)
# ---------------------------------------------------------------------
# Estimation of ATE 
# ---------------------------------------------------------------------
estimates <- tar_map(
    # outer map generates data under varying parameters of the data generating model
    values = varying,
    unlist = FALSE,
    tar_target(ESTIMATES,{
        ## lavaModel;
        out = do.call("rbind",mclapply(1:reps,function(b){
            set = FIXED
            set$formula.list = formula_settings[[setting]]
            simulated_data <- simulateData(setting = set,
                                           A1_T1 = A1_T1,
                                           A1_T2 = A1_T2,
                                           A2_T1 = A2_T1,
                                           A2_T2 = A2_T2,
                                           n = sample.size,
                                           scale.censored = scale.censored,
                                           keep.latent = FALSE)
            if (treat == "all"){
                ff = Hist(time,event)~intervene(A1)+intervene(A2)+intervene(A3)+intervene(A4)+intervene(A5)+intervene(A6)+intervene(A7)+intervene(A8)+intervene(A9)+intervene(A10)+X1+X2+X3+X4+X5+X6+X7
            } else{
                ff = Hist(time,event)~intervene(A1)+A2+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6+X7
            }
            if (weighter == "km"){
                x <- causalhunter(formula=ff,
                                  method = method,
                                  weighter=weighter,
                                  args.weight = list(num.trees = num.trees,alpha = 0.05,mtry = 17),
                                  fit.separate = TRUE,
                                  num.trees=num.trees,
                                  CR.as.censoring = net,
                                  data=simulated_data,
                                  times=horizon,
                                  formula.weight = Hist(time,event)~A1+A2)
            }else{
                x <- causalhunter(formula=ff,
                                  method = method,
                                  weighter=weighter,
                                  args.weight = list(num.trees = num.trees,alpha = 0.05,mtry = 17),
                                  fit.separate = TRUE,
                                  num.trees=num.trees,
                                  CR.as.censoring = net,
                                  data=simulated_data,
                                  times=horizon,
                                  formula.weight = Hist(time,event)~A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6+X7)
            }
            if (treat == "all"){
                x[,rank := rank(-abs(ate))]
            }else{
                x[,rank := NA]
            }
            x <- cbind(x,data.table(n = sample.size,
                                    net = as.numeric(net),
                                    scale.censored = scale.censored,
                                    num.trees = num.trees,
                                    method = method,
                                    A1_T1 = A1_T1,
                                    A1_T2 = A1_T2,
                                    A2_T1 = A2_T1,
                                    A2_T2 = A2_T2,
                                    horizon = horizon,
                                    formula = setting,
                                    weighter = weighter,
                                    theme = theme))
            x
        },mc.cores = MCCORES)) # end of repetitions
        gc()
        out
    },cue = tar_cue(mode = "never")) # end of parameter map
)

# combine
ate <- tar_combine(ESTIMATE_ATE,{
    estimates
})
# ---------------------------------------------------------------------
# Summarize performance of estimators against true parameter values
# ---------------------------------------------------------------------
results <- tar_target(RESULTS,
                      summarizePerformance(truth = TRUTH,
                                           estimate = ESTIMATE_ATE),
                      deployment = "main")


ranking <- tar_target(RANKING,
                      rankingPerformance(estimate = ESTIMATE_ATE),
                      deployment = "main")


plotframe <- tar_target(PLOTFRAME,{
    e <- ESTIMATE_ATE[intervene == "A1"]
    t <- TRUTH
    ## e <- tar_read("ESTIMATE_ATE")[intervene == "A1"]
    ## t <- tar_read("TRUTH")[intervene == "A1"]
    ## setnames(e,"time","horizon")
    setkeyv(e,c("net","formula","intervene","horizon","scale.censored","A1_T1","A1_T2","A2_T1","A2_T2"))
    setkeyv(t,c("net","formula","intervene","horizon","scale.censored","A1_T1","A1_T2","A2_T1","A2_T2"))
    t.ate = t[cause == 1 & intervene == "A1",.(true.ate = mean(ate)),keyby = c("net","formula","horizon","A1_T1","A1_T2","A2_T1","A2_T2")]
    t.cens = t[cause == 1 & intervene == "A1",.(censored.tau = mean(censored.tau)),keyby = c("formula","horizon","scale.censored")]
    e = merge(e,t.ate,by = c("net","formula","horizon","A1_T1","A1_T2","A2_T1","A2_T2"),all.y = FALSE)
    e = merge(e,t.cens,by = c("formula","horizon","scale.censored"),all.y = FALSE)
    e[,n:=factor(n)]
    e[,net:=factor(net)]
    e[,method:=factor(method,levels=c("causal_forest","CSC","FGR"),labels=c("Causal forest","CSC","FGR"))]
    e[,num.trees:=factor(num.trees)]
    e[,horizon:=factor(horizon)]
    e[,censored.tau:=factor(round(censored.tau,1))]
    e[,A1_T1:=factor(A1_T1)]
    e[,A1_T2:=factor(A1_T2)]
    e[,A2_T1:=factor(A2_T1)]
    e[,A2_T2:=factor(A2_T2)]
    e
})



######################################################################
### simulation_targets.R ends here
