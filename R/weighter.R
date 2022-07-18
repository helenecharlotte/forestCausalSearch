### weighter.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jun 27 2020 (06:33) 
## Version: 
## Last-Updated: Jun 23 2022 (09:39) 
##           By: Thomas Alexander Gerds
##     Update #: 198
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param formula
##' @param data
##' @param CR.as.censoring
##' @param cause
##' @param method
##' @param times
##' @param ...
##' @return 
##' @seealso 
##' @examples
##' library(riskRegression)
##' library(survival)
##' library(prodlim)
##' library(data.table)
##' library(grf)
##' library(ranger)
##' mydata <- sampleData(200,outcome="competing.risk")
##' weighter(formula=Hist(time,event)~X1+X5+X7,data=mydata,times=5)
##' set.seed(8)
##' weighter(formula=Hist(time,event)~X1+X5+X7,data=mydata,times=c(5,8),truncate=FALSE)
##' set.seed(8)
##' u=weighter(formula=Hist(time,event)~X1+X5+X7,data=mydata,times=c(5,8),truncate=TRUE)
##' v=weighter(formula=Hist(time,event)~X1+X5+X7,data=mydata,times=c(5,8),truncate=TRUE,method="km")
##' weighter(formula=Hist(time,event)~X1+X5+X7,data=mydata,times=c(5,8),truncate=FALSE,method="km")
##' a <- weighter(formula=Hist(time,event)~X1+X5+X7,data=mydata,times=5)
##' b <- weighter(formula=Hist(time,event)~X1+X5+X7,data=mydata,times=5,CR.as.censoring=TRUE)
##' cbind(a,b)
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
weighter <- function(formula,
                     data,
                     CR.as.censoring,
                     fit.separate=FALSE,
                     cause=1,
                     method="km",
                     times,
                     truncate=FALSE,
                     ...){
    # parse formula in data
    EHF <- prodlim::EventHistory.frame(formula=formula,data=data,specials=NULL,unspecialsDesign=FALSE)
    X <- EHF$design
    response.type <- attr(EHF$event.history,"model") # either "survival" or "competing.risk"
    if(response.type!="competing.risks") {
        if (CR.as.censoring[[1]] == TRUE){
            warning("No competing risks in data")
            CR.as.censoring <- FALSE
        }
        fit.separate <- FALSE
    }
    cens.type <- attr(EHF$event.history,"cens.type") # either "uncensored" or "rightCensored"
    dt <- data.table(cbind(unclass(EHF$event.history),X))
    #
    # in dt the response is time, status for survival
    # and time, event (1,2,3=censored), status (1,0=censored) for competing risks
    # in both cases we keep only the event variable (0=censored,1) or (0=censored,1,2)
    if (response.type=="survival") {
        setnames(dt,"status","event")
    } else{
        if (response.type=="competing.risks"){
            ## in competing risks case set censored to have event value 0
            dt[,event:=event*status]
            dt[,status := NULL]
        }
    }

    ## truncate all times that are larger than the maximum evaluation time
    if (truncate[[1]]==TRUE){
        if (any(time>max(times))){
            cens.type = "rightCensored"
            if (response.type=="survival")
                dt[time>max(times),event:=0]
            else
                dt[time>max(times),event:=0]
            dt[time>max(times),time:=max(times)+max(times)/10000]
        }
    }
    # case: uncensored and not intervening on competing risk
    if (CR.as.censoring==FALSE && cens.type=="uncensored"){
        Weight <- as.numeric(dt[["event"]] == cause)
    }else{
        ## cases: either right censored data or net effects in the hypothetical world
        ## where the competing risk has been removed
        dt[,Weight:=1]
        causes <- c(0,1,2)
        stopifnot(cause%in%causes)
        if (CR.as.censoring){
            if (fit.separate) {
                JJ <- causes[causes != cause]} else {JJ <- cause}
        }else{
            JJ <- cause
        }
        if (tolower(method)=="km"){
            for (jj in JJ){
                if (jj == cause){
                    if (CR.as.censoring){
                        # reverse Kaplan-Meier for combined censoring/cause2 distribution 
                        dt[,Status := as.numeric(event == cause)]
                    }else{
                        # reverse Kaplan-Meier for censoring distribution
                        dt[,Status := as.numeric(event)]}
                } else{
                    # separate reverse Kaplan-Meiers for censoring (jj=0)
                    # and latent cause 2 event time distribution (jj=2) 
                    dt[,Status := as.numeric(event != jj)]
                }
                ff <- update(formula,"Hist(time,Status)~.")
                reverse.km.fit <- prodlim(ff, data=dt, reverse=TRUE)
                ## result of predictSurvIndividual is sorted by (time, -status)
                dt[,id:=1:.N]
                setorder(dt,time,Status)
                dt[,Weight:=Weight*(event==cause)/predictSurvIndividual(reverse.km.fit)]
                ## back to original order
                setorder(dt,id)
            }
            Weight <- dt[["Weight"]]
        } else {
            # method: ranger
            for (jj in JJ) {
                if (jj == cause){
                    if (CR.as.censoring){
                        # forest combined censoring/cause2 distribution 
                        dt[,Status := as.numeric(event != cause)]
                    }else{
                        # forest for censoring distribution
                        dt[,Status := as.numeric(event == 0)]}
                } else{
                    # fit separate forests
                    # for censoring (jj=0)
                    # and latent cause 2 event time distribution (jj=2) 
                    dt[,Status := as.numeric(event == jj)]
                }
                ff <- update(formula,paste0("Surv(time,Status)~."))
                reverse.forest <- do.call("ranger",c(list(formula=ff,
                                                          data=dt[,all.vars(ff),with = FALSE],
                                                          replace = FALSE,
                                                          splitrule = "maxstat"),
                                                     ...))
                Gmat <- stats::predict(reverse.forest,data=dt)$survival
                jtimes <- ranger::timepoints(reverse.forest)
                Gi.minus <- sapply(1:length(dt$time),function(i){
                    pos <- prodlim::sindex(jump.times=jtimes,eval.times=dt$time[[i]],
                                           comp="smaller",strict=1L)
                    c(1,Gmat[i,])[1+pos]
                })
                rm(Gmat)
                gc()
                # multiply the weights G and G_2
                dt[,Weight:=Weight*(event == cause)/Gi.minus]
            }
            Weight <- dt[["Weight"]]
        }
    }
    Y <- do.call("cbind",lapply(times,function(t0){
        # weigth is already zero for event!=cause
        y = as.numeric(dt[["time"]]<=t0)*Weight
        if (any(is.infinite(y))) stop(paste0("Weighted outcome status has infinite values at time ",t0))
        if (length(unique(y))<=1) {
            print(table(y))
            print(head(dt[["time"]]))
            print(head(Weight))
            stop(paste0("Outcome status has no variation at time ",t0))
        }
        y
    }))
    Y
}


######################################################################
### weighter.R ends here
