### causalhunter.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jun  4 2020 (16:37) 
## Version: 
## Last-Updated: Jul 14 2022 (09:26) 
##           By: Thomas Alexander Gerds
##     Update #: 223
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
##' @param method One of \code{"causal_forests"}, \code{"FGR"), \code{"CSC"}, or \code{"naive"}
##' @param weighter Method for estimating censoring weights.
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
##' mydata <- sampleData(200,outcome="competing.risk")
##' cbind(weighter(Hist(time,event)~X1+X5+X7,data=mydata,times=5,method="km"),
##'       weighter(Hist(time,event)~X1+X5+X7,data=mydata,times=5,method="ranger"))
##' causalhunter(Hist(time,event)~intervene(X1)+X5+X7,data=mydata,times=5,weighter="ranger")
##' causalhunter(Hist(time,event)~intervene(X1)+X5+X7,data=mydata,times=5,method="naive")
##' causalhunter(Hist(time,event)~intervene(X1)+X5+X7,data=mydata,times=5,method="FGR")
##' causalhunter(Hist(time,event)~intervene(X1)+X5+X7,data=mydata,times=5,method="CSC")
##' # survival case
##' causalhunter(Hist(time,event!=0)~intervene(X1)+X5+X7,data=mydata,times=5,method="causal_forest")
##'
##' causalhunter(Hist(time,event)~intervene(X1)+intervene(X2)+X5+X7,data=mydata,times=5)
##' causalhunter(Hist(time,event)~intervene(X1)+intervene(X2)+X5+X7,data=mydata,times=c(5,8))
##' 
##' @export 
##' @author Helene Charlotte Rytgaard <hely@@sund.ku.dk>, Thomas A. Gerds <tag@@biostat.ku.dk>
causalhunter <- function(formula,
                         data,
                         CR.as.censoring=FALSE,
                         cause=1,
                         method = "causal_forest",
                         weighter="ranger",
                         fit.separate=FALSE,
                         formula.weight,
                         args.weight=NULL,
                         times,
                         truncate=FALSE,
                         ...){
    requireNamespace("data.table")
    requireNamespace("prodlim")
    Hist <- prodlim::Hist
    data <- data.table::copy(data)
    EHF <- prodlim::EventHistory.frame(formula=formula,
                                       data=data,
                                       stripSpecials="intervene",
                                       specials="intervene",
                                       specialsFactor=FALSE,
                                       unspecialsDesign=TRUE)
    # unspecialsDesign = TRUE because grf does not work with non-numeric features
    variables <- unique(c(c(names(attr(EHF$design, "levels")),
                            colnames(EHF$design))[c(names(attr(EHF$design, "levels")),
                                                    colnames(EHF$design))%in%colnames(data)],
                          colnames(EHF$intervene)))
    if (missing(formula.weight) || is.null(formula.weight)){
        # add all variables if formula.weight is missing
        formula.weight <- formula(paste0("Hist(time,event)~",paste(variables,collapse="+")))
    }
    # restrict to variables used
    dt <- data.table(cbind(unclass(EHF$event.history),
                           #EHF$design 
                           data[, variables, with=FALSE],
                           EHF$intervene))
    if (attr(EHF$event.history,"model")=="survival"){
        # no competing risks
        dt[,event:=status]
    } else{
        ## fix censored event status
        dt[status==0,event:=0]
    }
    out <- NULL
    switch(method,"causal_forest" = {
        Y <- do.call("weighter",
                     c(list(formula=formula.weight,data=dt,times=times,
                            method=weighter,CR.as.censoring=CR.as.censoring,
                            fit.separate=fit.separate,
                            cause=cause,truncate=truncate),args.weight))
        #-- run grf:
        for (s in 1:length(times)){
            for (i in 1:NCOL(EHF$intervene)){
                thisdrug = colnames(EHF$intervene)[i]
                # use all other treatments and covariates
                XX = cbind(EHF$design,EHF$intervene[,setdiff(colnames(EHF$intervene),thisdrug)])
                grf.A <- do.call(causal_forest,c(list(X=XX,
                                                      Y=Y[,s,drop=TRUE],
                                                      W=as.numeric(as.character(EHF$intervene[[i]]))),
                                                 ...))
                est <- average_treatment_effect(grf.A)
                names(est) <- c("ate","se")
                lower=est[[1]]+qnorm(0.025)*est[[2]]
                upper=est[[1]]+qnorm(0.975)*est[[2]]
                est <- data.frame(time=times[s],intervene=colnames(EHF$intervene)[[i]],ate=est[[1]],se=est[[2]],lower=lower,upper=upper)
                out <- rbind(out,est)
            }
        }
        rownames(out) <- NULL
        out
    },
    "naive" = {
        for (i in 1:NCOL(EHF$intervene)){
            A = names(EHF$intervene)[i]
            formula.i = update(formula,paste0(".~",A))
            fit.i = prodlim(formula.i,data = dt)
            nd = data.frame(c(0,1))
            names(nd) = A
            sum.i = summary(fit.i,newdata = nd,cause = cause,asMatrix = FALSE,times = times,type = "risk")$table[[cause]]
            ate.i = sum.i[[2]][,"cuminc"]-sum.i[[1]][,"cuminc"]
            se.i = sqrt((sum.i[[1]][,"se.cuminc"])^2+(sum.i[[2]][,"se.cuminc"])^2)
            lower.i = ate.i+qnorm(0.025)*se.i
            upper.i = ate.i+qnorm(0.975)*se.i
            est = data.table::data.table(time = times,intervene = A,ate = ate.i,se = se.i,lower = lower.i,upper = upper.i)
            out <- rbind(out,est)
        }
    },
    "FGR" = {
        stopifnot(attr(EHF$event.history,"model") == "competing.risks")
        for (i in 1:NCOL(EHF$intervene)){
            A = names(EHF$intervene)[i]
            formula.i = update(formula,paste0(".~",paste(attr(attr(EHF,"Terms"),"term.labels"),collapse = "+")))
            fit.i <- FGR(formula.i,data = data,cause = cause)
            data0 <- data1 <- copy(dt)
            data0[[A]] <- 0
            data1[[A]] <- 1
            ate.i <- mean(predictRisk(fit.i,newdata=data1,times=times)- predictRisk(fit.i,newdata=data0,times=times))
            est = data.table::data.table(time = times,intervene = A,ate = ate.i,se = NA,lower = NA,upper = NA)
            out <- rbind(out,est)
        }
    },"CSC" = {
        stopifnot(attr(EHF$event.history,"model") == "competing.risks")
        for (i in 1:NCOL(EHF$intervene)){
            AA = names(EHF$intervene)[i]
            formula.i = update(formula,paste0(".~",paste(attr(attr(EHF,"Terms"),"term.labels"),collapse = "+")))
            data[[AA]] = factor(data[[AA]])
            fit.i = CSC(formula.i,data = data,cause = cause)
            fit.i$call$formula = eval(formula.i)
            res = ate(fit.i,data = data,treatment = AA,times = times,verbose = FALSE)$diffRisk[,data.table(time = time,intervene = AA,ate = estimate, se = se,lower = lower,upper = upper)]
            out <- rbind(out,res)
        }
    })
    data.table::setDT(out)
    out
}

######################################################################
### causalhunter.R ends here
