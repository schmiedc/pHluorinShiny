library(tidyverse)
library(broom)
library(plyr)

calcTau <- function(inputDataframe) {
    
        # inputDataframe = finalTable
        
        count <- as.data.frame(table(inputDataframe$name))
        
        tau.list <- list()
        
        for (names in count$Var1){
            
            name.table <-subset(inputDataframe, name == names )
            
            # detect peak of curve, get time and then cut off data at that point
            reuptake <- subset( name.table, time >= name.table$time[ which.max( name.table$peak_norm ) ] )
            
            # self-starting function, 
            # a special function for curve fitting that guesses its own start parameters. 
            # The asymptotic regression function, SSasymp is equivalent to our exponential decay:
            
            tau <- "NaN"
            rse <- "NoFit"
            degreeFreedom <- "NaN"
            
            tryCatch({
            
                fit <- nls(peak_norm ~ SSasymp(time, yf, y0, log_alpha), data = reuptake)
                
                #fit <- nls(peak_norm ~ yf + (y0 - yf) * exp(-alpha * time), data = reuptake, start = list(y0 = 0, yf = 0, alpha = 0))
                
                
                summary <- summary(fit)
                rse <- summary$sigma
                df <- summary$df
                degreeFreedom <- df[2]
            
                # summary(fit)
                fittingParams <- coef(fit)
            
                yf = fittingParams[1]
                y0 = fittingParams[2]
                alpha = exp(fittingParams[3])
                y = 1 / exp(1)
            
                # solve for t (timeconstant) if y is given
                tau = log( (y - yf) / (y0 - yf) ) / - alpha
            
                #tauDataframe <- data.frame(".id" = reuptake$.id[1], "name" = reuptake$name[1], "tau" = tau)
                
            }, warning = function(w) {
                print("WARNING: Fitting warning produced")
            }, error = function(e) {
                
                print("ERROR: No fit produced")
                
                tau <- "NaN"
                rse <- "NoFit"
                degreeFreedom <- "NaN"
                
            }, finally = {
                
                tau.list[[names]]$tau <- tau
                tau.list[[names]]$rse <- rse
                tau.list[[names]]$df <- degreeFreedom
                
            })
                 
        }
        
        tau <- ldply(tau.list, data.frame)
        names(tau)[1] <- "name"
        return (tau)
        gc()

}
    