library(tidyverse)
library(broom)

dir = "/data1/FMP_Docs/Repositories/minimal-datasets_FMP/pHlorin_TS/ROutput/"
curve <- read.csv(paste0(dir ,"Mean_Test.csv"), header =TRUE)

# detect peak of curve, get time and then cut off data at that point
reuptake <- subset( curve, time >= curve$time[ which.max( curve$peak_norm ) ] )
plot(reuptake$time, reuptake$peak_norm)

head(reuptake)

# self-starting function, 
# a special function for curve fitting that guesses its own start parameters. 
# The asymptotic regression function, SSasymp is equivalent to our exponential decay:
fit <- nls(peak_norm ~ SSasymp(time, yf, y0, log_alpha), data = reuptake)

summary(fit)
summary <- summary(fit)
summary$sigma
summary$df[2]

coef(fit)
fittingParams <- coef(fit)
lines(reuptake$time,predict(fit))

?summary.nls

alpha = exp(fittingParams[3])
alpha
yf = fittingParams[1]
y0 = fittingParams[2]
y = 1 / exp(1)

yf
y0

# solve for t (timeconstant) if y (yTimeConstant) is given
# yTimeConstant = yf + (y0 - yf) * exp(-alpha * timeconstant)
# yTimeConstant - yf = (y0 - yf) * exp(-alpha * timeconstant)
# ( yTimeConstant - yf ) / (y0 - yf) = exp(-alpha * timeconstant)
# log( ( yTimeConstant - yf ) / (y0 - yf) ) = -alpha * timeconstant
tau = log( (y - yf) / (y0 - yf) ) / - alpha
tau

qplot(reuptake$time, reuptake$peak_norm, data =  augment(fit)) + 
    geom_hline(yintercept=y) +
    geom_vline(xintercept=tau) +
    geom_line(aes(y = .fitted) 
    )

# fitting with nls
estFit <- nls(
    #  function we use to fi
    peak_norm~ yf + (y0 - yf) * exp(-alpha * time), 
    # inputdata
    data = reuptake,
    # starting values of fit
    start = list(
        y0 = 3.789088,
        yf = 0.08959001, 
        alpha = alpha))


summary(estFit)
t = 26.88
yTimeConstant = 0.3679
y = yf + (y0 - yf) * exp(-alpha * t)

