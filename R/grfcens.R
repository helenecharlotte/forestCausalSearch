### grfcens.R --- 
#----------------------------------------------------------------------
## Author: Helene Charlotte Wiese Rytgaard & Thomas Alexander Gerds
## Created: Jun  4 2020 (16:37) 
## Version: 
## Last-Updated: May 11 2022 (06:59) 
##           By: Thomas Alexander Gerds
##     Update #: 107
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' This function implements a two-step estimator of the average treatment effect 
##' for time-to-event outcome defined as a t-year risk difference
##' according to an intervention. In the presence of competing risks, the function
##' estimates either net or crude probabilities.
##' 
##' One should check the causal three: consistency, positivity and exchangeability.
##' Exchangeability may be a more serious concern when analysing the net probabilities.
##' @title Causal forests for censored data
##' @param formula Formula where left hand side specifies the event history variables
##' and the right hand side the predictor variables where exactly one variable is
##' marked as the exposure of interest via intervene(variable). See examples.
##' @param data Data in which to evaluate the formulae.
##' @param CR.as.censoring Only for settings with competing risks. Logical.
##' If \code{TRUE} analyse net probabilities.
##' @param cause Only for settings with competing risks. The cause of interest.
##' @param method.weight Method for estimating censoring weights.
##' Either \code{'km'} or \code{'ranger'}.
##' @param formula.weight Optional formula specifying the variables used to
##' estimate the censoring weights.
##' @param args.weight Additional arguments passed to method.weight
##' for estimating the censoring weigts.
##' @param times Time point(s) for evaluating the t-year risks.
##' @param ... Additional arguments passed to \code{causal_forest}.
##' @return Average treatment effects
##' @seealso causal_forest
##' @examples
##' library(lava)
##' library(riskRegression)
##' library(prodlim)
##' library(riskRegression)
##' library(prodlim)
##' library(survival)
##' library(data.table)
##' library(grf)
##' library(ranger)
##' m <- lvm(~A+X1+X2+X3+X4+X5)
##' # treatment variable
##' distribution(m,~A) <- binomial.lvm()
##' # latent outcome variables
##' distribution(m,~latenttime1) <- coxWeibull.lvm()
##' distribution(m,~latenttime2) <- coxWeibull.lvm()
##' distribution(m,~censtime) <- coxWeibull.lvm()
##' # observed outcome
##' m <- eventTime(m,time~min(latenttime1=1,latenttime2=2,censtime=0),"event")
##' # dependencies: log-odds ratios and log-hazard ratios
##' # no effect of A on latenttime1
##' regression(m,A~X1+X2+X3+X4+X5) <- c(.1,-.3,.8,-.1,0)
##' regression(m,latenttime1~A+X1+X2+X3+X4+X5) <- c(0,.1,-.3,.8,-.1,0)
##' regression(m,latenttime2~A+X1+X2+X3+X4+X5) <- c(0.2,0,-.3,-.8,-.1,.5)
##' # no effect of any variable on censoring distribution
##' gammaA <- gamma1 <- gamma2 <- gamma3 <- gamma4 <- gamma5  <- 0
##' regression(m,censtime~A+X1+X2+X3+X4+X5) <- c(gammaA,gamma1,gamma2,gamma3,gamma4,gamma5)
##' plot(m)
##' set.seed(99)
##' d <- sim(m,200)
##' plot(prodlim(Hist(time,event)~A,data=d))
##' fit0 <- prodlim(Hist(time,event)~A,data=d)
##' fit1 <- grfcens(formula=Hist(time,event)~intervene(A)+X1+X2+X3+X4+X5,data=d,times=7)
##' fit1
##' fit2 <- grfcens(formula=Hist(time,event)~intervene(A)+X1+X2+X3+X4+X5,data=d,times=7,
##'                 num.tree=50,args.weight=list(num.tree=50))
##' fit2
##' # with effect of A on T1
##' regression(m,latenttime1~A) <- 0.5
##' set.seed(98)
##' d <- sim(m,200)
##' fit2a <- grfcens(formula=Hist(time,event)~intervene(A)+X1+X2+X3+X4+X5,data=d,times=7,
##'                 num.tree=50,args.weight=list(num.tree=50))
##' fit2a
##' # with effect of A on T2 in opposite direction
##' regression(m,latenttime2~A) <- -0.9
##' set.seed(98)
##' d <- sim(m,200)
##' fit2b <- grfcens(formula=Hist(time,event)~intervene(A)+X1+X2+X3+X4+X5,data=d,times=7,
##'                 num.tree=50,args.weight=list(num.tree=50))
##' fit2b
##' 
##' ## CR.as.censoring:
##' # a) variable A has an effect only through the rate of the competing risk
##' regression(m,latenttime1~A) <- 0
##' regression(m,latenttime2~A) <- 0
##' set.seed(89)
##' x <- foreach(i=1:10,.combine="rbind")%dopar%{
##' d <- sim(m,200)
##' fit3  <- grfcens(formula=Hist(time,event)~intervene(A)+X1+X2+X3+X4+X5,data=d,times=7,
##'                  num.tree=50,args.weight=list(num.tree=50),CR.as.censoring=FALSE)
##' fit3b  <- grfcens(formula=Hist(time,event)~intervene(A)+X1+X2+X3+X4+X5,data=d,times=7,
##'                   num.tree=50,args.weight=list(num.tree=50),CR.as.censoring=TRUE)
##' c(here=fit3[[1]],hypo=fit3b[[2]])
##' }
##' 
##' @export 
##' @author Helene Charlotte Wiese Rytgaard <hely@@sund.ku.dk>, Thomas A. Gerds <tag@@biostat.ku.dk>
grfcens <- function(formula,
                    data, 
                    CR.as.censoring=FALSE,
                    cause=1,
                    method.weight="ranger",
                    formula.weight,
                    args.weight=NULL,
                    times,...){
    EHF <- prodlim::EventHistory.frame(formula=formula,
                                       data=data,
                                       stripSpecials="intervene",
                                       specials="intervene")
    if (length(times)>1) warning("Works currently only for one time point at a time.")
    t0=times[[1]]
    response <- EHF$event.history
    response.type <- attr(response,"model") # either "survival" or "competing.risk"
    cens.type <- attr(EHF$event.history,"cens.type") # either "uncensored" or "rightCensored"
    if (CR.as.censoring == TRUE && response.type=="survival")
        warning("Only one cause of event in data. Removing a competing risk which does not occur.")
    # uncensored
    if (CR.as.censoring==FALSE &&  cens.type=="uncensored"){
        Y <- as.numeric(response[,"time"]<=t0)
    }else{
        ## require either right censored data, or request for the hypothetical world
        ## where competing risk has been removed
        stopifnot(cens.type=="rightCensored" || CR.as.censoring == TRUE)
        dt <- data.table(cbind(unclass(EHF$event.history),EHF$design,EHF$intervene))
        ## in case without competing risks add event variable
        if (response.type=="survival") dt[,event:=status]
        if (CR.as.censoring) {
            dt[,is.censored:=as.numeric(event!=cause)]
        }else{
            dt[,is.censored:=as.numeric(status==0)]
        }
        if (tolower(method.weight)=="km"){
            if (CR.as.censoring==TRUE) dt[is.censored==1,status:=0]
            if (missing(formula.weight)){
                formula.weight <- formula(paste0("Hist(time,status)~",paste(c(colnames(EHF$design),colnames(EHF$intervene)),collapse="+")))
                reverse.km.fit <- prodlim(formula.weight, data=dt, reverse=TRUE)
            }else{
                formula.weight <- update.formula(formula.weight,"Hist(time,status)~.")
                reverse.km.fit <- prodlim(formula.weight, data=dt, reverse=TRUE)
            }
            ## result of predictSurvIndividual is sorted by (time, -status)
            dt[,id:=1:.N]
            setorder(dt,time,is.censored)
            #-- define real-valued outcome:
            dt[,Y:=(event==cause)*(time<=t0)/predictSurvIndividual(reverse.km.fit)]
            ## back to original order
            setorder(dt,id)
            Y <- dt[["Y"]]
        }else{
            if (missing(formula.weight)){
                ff <- formula(paste0("Surv(time,is.censored)~",paste(c(colnames(EHF$design),colnames(EHF$intervene)),collapse="+")))
            }else{
                ff <- update(formula.weight,"Surv(time,is.censored)~.")
            }
            reverse.forest <- do.call("ranger",c(list(formula=ff, data=dt),args.weight))
            Gmat <- stats::predict(reverse.forest,data=dt)$survival
            jtimes <- ranger::timepoints(reverse.forest)
            Gi.minus <- sapply(1:length(dt$time),function(i){
                pos <- prodlim::sindex(jump.times=jtimes,
                                       eval.times=dt$time[[i]],
                                       comp="smaller",
                                       strict=1L)
                c(1,Gmat[i,])[1+pos]
            })
            dt[,weights:=Gi.minus]
            dt[,Y:=as.numeric((event==cause)*(time<=t0))]
            dt[Y!=0,Y:=1/weights]
            Y <- dt[["Y"]]
        }
    }
    if (any(is.infinite(Y))) stop(paste0("Weighted outcome status has infinite values at time ",t0))
    if (length(unique(Y))<=1) stop(paste0("Outcome status has no variation at time ",t0))
    #-- run grf:
    grf.A <- causal_forest(X=EHF$design,
                           Y=Y,
                           W=as.numeric(as.character(EHF$intervene[[1]])),
                           ...)
    out <- average_treatment_effect(grf.A)
    names(out) <- c("ate","se")
    lower=out[[1]]+qnorm(0.025)*out[[2]]
    upper=out[[1]]+qnorm(0.975)*out[[2]]
    out <- c(out,lower=lower,upper=upper)
    out
}

######################################################################
### grfcens.R ends here
