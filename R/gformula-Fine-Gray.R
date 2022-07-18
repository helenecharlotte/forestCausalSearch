### gformula-Fine-Gray.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jan 14 2022 (09:26) 
## Version: 
## Last-Updated: Jan 14 2022 (10:10) 
##           By: Thomas Alexander Gerds
##     Update #: 1
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
gformulaFG <- function(fit,data,intervene,times){
    data0 <- data1 <- copy(data)
    data0[[intervene]] <- 0
    data1[[intervene]] <- 1
    ate <- mean(predictRisk(fit,newdata=data1,times=times)- predictRisk(fit,newdata=data0,times=times))
    ate
}


#----------------------------------------------------------------------
### gformula-Fine-Gray.R ends here
