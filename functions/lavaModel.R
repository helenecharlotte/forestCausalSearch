### lavaModel.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr  9 2022 (19:27) 
## Version: 
## Last-Updated: Jul 11 2022 (10:23) 
##           By: Thomas Alexander Gerds
##     Update #: 65
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
lavaModel <- function(event.times = c("T1","T2"),
                      treatments,
                      binary.covariates = NULL,
                      normal.covariates = NULL,
                      quadratic.covariates = NULL,
                      scale = 1/100,
                      treatment.effect,
                      formula.list = NULL,
                      randomized = FALSE,
                      ...){
    # empty latent variable model
    requireNamespace("lava")
    m <- lava::lvm()
    # covariates
    for (X in names(binary.covariates)) {
        lava::distribution(m,X) <- lava::binomial.lvm()
        lava::intercept(m,X) <- log(binary.covariates[[X]])
    }
    for (X in normal.covariates) lava::distribution(m,X) <- lava::normal.lvm()
    for (X in quadratic.covariates) lava::transform(m,formula(paste0(X,"_2","~",X))) <- function(x){x^2}
    # treatment variables
    for (A in names(treatments)) {
        lava::distribution(m,A) <- lava::binomial.lvm()
        lava::intercept(m,A) <- log(treatments[[A]])
    }
    # initialize all latent event times with same scale
    latent.times = sapply(event.times,function(x){paste0(x,"_",c("placebo_placebo","treated_placebo","placebo_treated","treated_treated"))})
    for (R in c(latent.times,"C")) lava::distribution(m,R) <- lava::coxWeibull.lvm(scale = scale)
    # change the scale according to treatment effects
    lava::distribution(m,"T1_treated_treated") <- lava::coxWeibull.lvm(scale = scale * treatment.effect[["A1_T1"]]* treatment.effect[["A2_T1"]])
    lava::distribution(m,"T1_treated_placebo") <- lava::coxWeibull.lvm(scale = scale * treatment.effect[["A1_T1"]])
    lava::distribution(m,"T1_placebo_treated") <- lava::coxWeibull.lvm(scale = scale * treatment.effect[["A2_T1"]])
    lava::distribution(m,"T2_treated_treated") <- lava::coxWeibull.lvm(scale = scale * treatment.effect[["A1_T2"]]* treatment.effect[["A2_T2"]])
    lava::distribution(m,"T2_treated_placebo") <- lava::coxWeibull.lvm(scale = scale * treatment.effect[["A1_T2"]])
    lava::distribution(m,"T2_placebo_treated") <- lava::coxWeibull.lvm(scale = scale * treatment.effect[["A2_T2"]])
    # all other regressions
    # all other regressions
    for (f in formula.list) {
        if ((timevar <- all.vars(f)[[1]]) %in% c("T1","T2")){
            for(a in c("_treated_treated","_treated_placebo","_placebo_treated","_placebo_placebo")){
                ## print(update(f,paste(paste0(timevar,a),"~.")))
                lava::regression(m) = update(f,paste(paste0(timevar,a),"~."))
            }
        }else{
            lava::regression(m) = f
        }
    }
    # event times under observed treatments
    lava::transform(m,T1~T1_placebo_placebo+T1_treated_placebo+T1_placebo_treated+T1_treated_treated+A1+A2) = function(x){
        x[["T1_placebo_placebo"]]*(x[["A1"]] == 0)*(x[["A2"]] == 0) +
            x[["T1_treated_placebo"]]*(x[["A1"]] == 1)*(x[["A2"]] == 0) +
            x[["T1_placebo_treated"]]*(x[["A1"]] == 0)*(x[["A2"]] == 1) +
            x[["T1_treated_treated"]]*(x[["A1"]] == 1)*(x[["A2"]] == 1)
    }
    lava::transform(m,T2~T2_placebo_placebo+T2_treated_placebo+T2_placebo_treated+T2_treated_treated+A1+A2) = function(x){
        x[["T2_placebo_placebo"]]*(x[["A1"]] == 0)*(x[["A2"]] == 0) +
            x[["T2_treated_placebo"]]*(x[["A1"]] == 1)*(x[["A2"]] == 0) +
            x[["T2_placebo_treated"]]*(x[["A1"]] == 0)*(x[["A2"]] == 1) +
            x[["T2_treated_treated"]]*(x[["A1"]] == 1)*(x[["A2"]] == 1)
    }
    # event times in hypothetical world where one of the treatments (A1, A2)
    # is as observed but the other one is set to a value 
    lava::transform(m,T1_placebo_A1~T1_placebo_placebo+T1_treated_placebo+T1_placebo_treated+T1_treated_treated+A2) = function(x){
        x[["T1_placebo_placebo"]]*(x[["A2"]] == 0) + x[["T1_placebo_treated"]]*(x[["A2"]] == 1) 
    }
    lava::transform(m,T1_treated_A1~T1_placebo_placebo+T1_treated_placebo+T1_placebo_treated+T1_treated_treated+A2) = function(x){
        x[["T1_treated_placebo"]]*(x[["A2"]] == 0) + x[["T1_treated_treated"]]*(x[["A2"]] == 1) 
    }
    lava::transform(m,T1_placebo_A2~T1_placebo_placebo+T1_treated_placebo+T1_placebo_treated+T1_treated_treated+A1) = function(x){
        x[["T1_placebo_placebo"]]*(x[["A1"]] == 0) + x[["T1_treated_placebo"]]*(x[["A1"]] == 1) 
    }
    lava::transform(m,T1_treated_A2~T1_placebo_placebo+T1_treated_placebo+T1_placebo_treated+T1_treated_treated+A1) = function(x){
        x[["T1_placebo_treated"]]*(x[["A1"]] == 0) + x[["T1_treated_treated"]]*(x[["A1"]] == 1) 
    }
    lava::transform(m,T2_placebo_A1~T2_placebo_placebo+T2_treated_placebo+T2_placebo_treated+T2_treated_treated+A2) = function(x){
        x[["T2_placebo_placebo"]]*(x[["A2"]] == 0) + x[["T2_placebo_treated"]]*(x[["A2"]] == 1) 
    }
    lava::transform(m,T2_treated_A1~T2_placebo_placebo+T2_treated_placebo+T2_placebo_treated+T2_treated_treated+A2) = function(x){
        x[["T2_treated_placebo"]]*(x[["A2"]] == 0) + x[["T2_treated_treated"]]*(x[["A2"]] == 1) 
    }
    lava::transform(m,T2_placebo_A2~T2_placebo_placebo+T2_treated_placebo+T2_placebo_treated+T2_treated_treated+A1) = function(x){
        x[["T2_placebo_placebo"]]*(x[["A1"]] == 0) + x[["T2_treated_placebo"]]*(x[["A1"]] == 1) 
    }
    lava::transform(m,T2_treated_A2~T2_placebo_placebo+T2_treated_placebo+T2_placebo_treated+T2_treated_treated+A1) = function(x){
        x[["T2_placebo_treated"]]*(x[["A1"]] == 0) + x[["T2_treated_treated"]]*(x[["A1"]] == 1) 
    }    
    # event time outcome
    m <- lava::eventTime(m,time~min(T1=1,T2=2,C=0),"event")
    m
}
## lava::distribution(m,"A1_0") <- lava::constant.lvm(value=0)
## lava::distribution(m,"A1_1") <- lava::constant.lvm(value=1)
## u <- lavaModel(treatments = c("A1","A2"),
               ## binary.covariates = paste0("X",1:3),
               ## normal.covariates = "Z",
               ## treatment.effect =list("A1_T1" = log(2), "A1_T2" = 0, "A2_T1" = 0, "A2_T2" = log(.5)),
               ## formula.list = NULL);sim(u,4)

## v <- lavaModel(seed = 11,horizon = 5,event.times = c("T1","T2"),
               ## treatments = paste0("A",1:10),
               ## binary.covariates = paste0("X",1:5),
               ## normal.covariates = paste0("X",6:7),
               ## treatment.effect = list("A1_T1" = log(2), "A1_T2" = 0,
                                       ## "A2_T1" = 0, "A2_T2" = log(.5)),
               ## formula.list =  list (T1 ~ f(X1,-.4) + f(X2,.3) + f(X6,-.1),
                                     ## T2  ~ f(X1,-.1) + f(X2,.6) + f(X6,.1),
                                     ## C ~ f(A1,.5)+ f(A3,.5) + f(X2,.6) + f(A6,.1),
                                     ## A1 ~ f(X1,.4) + f(X6,.8)),n = 4000)
## sim(v,4)
######################################################################
### lavaModel.R ends here
