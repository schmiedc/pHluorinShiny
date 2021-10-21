setwd("/data1/FMP_Docs/Repositories/plugins_FMP/pHluorinShiny/")
library(gridExtra)
source("dataProcessing.R")
source("saveData.R")
source("plotData.R")
source("fitting.R")

# ============================================================================
#
#  DESCRIPTION: Analyse Manual vs Automatic from processed results
#               Ctrl vs Treatment
#              
#       AUTHOR: Christopher Schmied, 
#      CONTACT: schmied@dzne.de
#     INSITUTE: Leibniz-Forschungsinstitut f r Molekulare Pharmakologie (FMP)
#               Cellular Imaging - Core facility
#               Campus Berlin-Buch
#               Robert-Roessle-Str. 10
#               13125 Berlin, Germany
#
#         BUGS:
#        NOTES: 
# DEPENDENCIES: ggplot2: install.packages("ggplot2")
#               xlsx: install.packages("gxlsx")
#               reshape2: install.packages("reshape2")
#               plyr: install.packages("plyr")
#               gridExtra: install.packages("gridExtra")
#               tidyverse: install.packages("tidyverse")
#               broom: install.packages("broom")
#
#      VERSION: 1.0.0
#      CREATED: 2018-05-24
#     REVISION: 2020-02-07
#
# ============================================================================
# where to get the files
indir = "/data1/FMP_Docs/Projects/Publication_SynapseJ/pHluorinJ_Data/RevisedAnalysis/"

# where to save the data
outdir = indir

# ============================================================================
resultname = "Test"

# Time resolution in seconds
timeResolution = 2

# when stimulation happens
# these frames are used for calcuating the mean intensity
# then this value is used for the surface normalization
frameStimulation = 5

# further settings
labelSignal = "Spot"
labelBackground = "background"

# get raw data
table.signal <- read_csv(paste0(indir,"_RawSignal.csv"))
table.background <- read_csv(paste0(indir,"_RawBackground.csv"))

# extracting experimental information from file name
table.signal <- table.signal %>% separate(name, 
                                          sep ="_", c("day", "treatment", "number"), 
                                          remove=FALSE)
table.background <- table.background %>% separate(name, 
                                                  sep ="_", c("day", "treatment", "number"), 
                                                  remove=FALSE)

# reduce data for number of ROIs and area
singleData_area <- subset(table.signal, variable == "area")
singleData_area <- subset(singleData_area, time == 0)

# compute & plot number of ROIs
roiNumber <- singleData_area %>% group_by(treatment) %>% summarize(count = n())

ggplot(data=roiNumber, aes(x=treatment, y=count)) +
  geom_bar(stat="identity") +
  ylab("Count") + 
  scale_y_continuous(limits = c(0, 6000), breaks = seq(0, 6000, by = 500)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

# compute average area of ROI
ggplot(data=singleData_area, aes(x=treatment, y=value)) +
  geom_boxplot(outlier.colour="black", outlitreatmenter.shape=16, outlier.size=2, notch=FALSE) +
  stat_boxplot(geom = 'errorbar', width = 0.1) +
  expand_limits(y = 0) +
  # scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, by = 1)) +
  ylab("Area (Micron)") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

# filter traces where peak is in the stimulation range (20s)
table.signal_mean_filter <- subset(table.signal, variable == "mean")
peaks <- table.signal_mean_filter %>% group_by(name, roi) %>% summarise(value = max(value))
peaks_frame <- left_join(peaks, table.signal_mean_filter, by = c("name", "roi", "value"))

filtered_peaks <- peaks_frame %>% filter(time <= 20)
filtered_peaks_2 <- filtered_peaks %>% select(c(-value, -frame, -time, -variable, -value) )

filtered_signal <- left_join(filtered_peaks_2, table.signal_mean_filter, by = c("name", "roi", "day", "treatment", "number"))

# calculates average mean intensity per frame per experiment
table.signal_mean <- subset(filtered_signal, variable == "mean")
table.signal_avg <- table.signal_mean %>% group_by(day, treatment, frame, time) %>% summarize(mean=mean(value), N = length(value), sd = sd(value), se = sd / sqrt(N))

table.background_mean <- subset(table.background, variable == "mean")
table.background_avg <- table.background_mean %>% group_by(day, treatment, frame, time) %>% summarize(mean=mean(value), N = length(value), sd = sd(value), se = sd / sqrt(N))

# generate final table
forBackgroundSubtraction <- merge(table.signal_avg, table.background_avg, by=c("day", "treatment", "frame", "time"), suffixes=c(".sig",".back"))

# normalize mean signal with mean background intensity
forBackgroundSubtraction$mean.corr <- forBackgroundSubtraction$mean.sig - forBackgroundSubtraction$mean.back

# ==============================================================================
# normalizations
forBackgroundSubtraction$name <- paste0(forBackgroundSubtraction$day, "_", forBackgroundSubtraction$treatment)

surfaceNormalized <- surfaceNormalisation(forBackgroundSubtraction, frameStimulation)

peakNormalized <- peakNormalisation(surfaceNormalized)

finalTable <- sortFrames(peakNormalized)

# ==============================================================================
# plots treatment vs ctrl
ggplot(data=finalTable, aes(x=time, y=surf_norm, group = name, color = name)) +
  geom_line() + 
  theme_light() +
  expand_limits(x = 0, y = 0.9) +
  # scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  # scale_y_continuous(limits = c(0.9, 2.5), breaks = seq(0.9, 2.5, by = 0.1)) +
  xlab("Time (s)") + 
  ylab("Norm. fluorescence intensity (A.U.)") + 
  ggtitle("Surf Norm ") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

ggplot(finalTable, aes(x=time, y=peak_norm, group = name, color = name)) +
  geom_line() + 
  theme_light() +
  # expand_limits(x = 0, y = 0.9) +
  # scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  # scale_y_continuous(limits = c(0.9, 1.7), breaks = seq(0.9, 1.7, by = 0.1)) +
  xlab("Time (s)") + 
  ylab("Norm. fluorescence intensity (A.U.)") + 
  ggtitle("Peak Norm") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

# ==============================================================================
finalTable_avg_surf <- finalTable %>% group_by(treatment, frame, time) %>% summarize(mean=mean(surf_norm), N = length(surf_norm), sd = sd(surf_norm), se = sd / sqrt(N))

finalTable_avg_surf$high_corr <- with(finalTable_avg_surf, finalTable_avg_surf$mean + finalTable_avg_surf$se)
finalTable_avg_surf$low_corr <-  with(finalTable_avg_surf, finalTable_avg_surf$mean - finalTable_avg_surf$se)

head(finalTable_avg_surf)

ggplot(finalTable_avg_surf, aes(x=time, y=mean, group = treatment, color = treatment)) +
  geom_line() + 
  guides(colour=FALSE)  + 
  geom_ribbon(aes(ymin = low_corr, ymax = high_corr, colour=treatment, group=treatment, fill = treatment ), alpha=.3) +
  # expand_limits(x = 0, y = 0.9) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(limits = c(0.9, 2.6), breaks = seq(0.9, 2.6, by = 0.1)) +
  xlab("Time (s)") + 
  ylab("Norm. fluorescence intensity (A.U.)") + 
  ggtitle("Avg. Surf Norm") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

finalTable_avg_peak <- finalTable %>% group_by(treatment, frame, time) %>% summarize(mean=mean(peak_norm), N = length(peak_norm), sd = sd(peak_norm), se = sd / sqrt(N))

finalTable_avg_peak$high_corr <- with(finalTable_avg_peak, finalTable_avg_peak$mean + finalTable_avg_peak$se)
finalTable_avg_peak$low_corr <-  with(finalTable_avg_peak, finalTable_avg_peak$mean - finalTable_avg_peak$se)

ggplot(finalTable_avg_peak, aes(x=time, y=mean, group = treatment, color = treatment)) +
  geom_line() + 
  guides(colour=FALSE)  + 
  geom_ribbon(aes(ymin = low_corr, ymax = high_corr, colour=treatment, group=treatment, fill = treatment ), alpha=.3) +
  # expand_limits(x = 0, y = 0.9) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  #scale_y_continuous(limits = c(0.8, 2.6), breaks = seq(0.8, 2.6, by = 0.2)) +
  xlab("Time (s)") + 
  ylab("Norm. fluorescence intensity (A.U.)") + 
  ggtitle("Avg. Peak Norm") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

# ==============================================================================
# compute peaks
peaks <- finalTable %>% group_by(name) %>% summarise(max = max(surf_norm))
peaks$deltaMax <- peaks$max - 1
peaks <- peaks %>% separate(name, sep ="_", c("day", "treatment"), remove=FALSE)

res_peaks <- t.test(deltaMax ~ treatment, data = peaks, paired = TRUE)
peaks

auto_peak_sd <- peaks %>% group_by(treatment) %>% summarize(mean=mean(deltaMax), N = length(deltaMax), sd = sd(deltaMax), se = sd / sqrt(N))
auto_peak_sd


# plot peak difference
ggplot(data=peaks, aes(x=treatment, y=deltaMax)) +
  geom_boxplot(outlier.size = 0, outlier.shape = 1) +
  stat_boxplot(geom = 'errorbar', width = 0.2) +
  geom_jitter(width = 0.1) +
  ylab("delta F (exocytosis)") + 
  expand_limits(y = 0) +
  scale_y_continuous(limits = c(0, 1.8), breaks = seq(0, 1.8, by = 0.2)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

# compute taus
tau <- calcTau(finalTable, 20)
tau <- tau %>% separate(name, sep ="_", c("day", "treatment"), remove=FALSE)

res_tau <- t.test(tau ~ treatment, data = tau, paired = TRUE)
res_tau

ggplot(data=tau, aes(x=treatment, y=tau)) +
  geom_boxplot(outlier.size = 0, outlier.shape = 1) +
  stat_boxplot(geom = 'errorbar', width = 0.2) +
  geom_jitter(width = 0.1) +
  ylab("tau") + 
  #expand_limits(y = 25) +
  #scale_y_continuous(limits = c(25, 125, breaks = seq(25, 125, by = 10)))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

# ==============================================================================
# manual analysis
table.signal_manual <- read_csv(paste0(indir,"ManualAnalysisResults.csv"))

# ==============================================================================
finalTable_avg_surf_manual <- table.signal_manual %>% group_by(treatment, time) %>% summarize(mean=mean(surf_norm), N = length(surf_norm), sd = sd(surf_norm), se = sd / sqrt(N))

finalTable_avg_surf_manual$high_corr <- with(finalTable_avg_surf_manual, finalTable_avg_surf_manual$mean + finalTable_avg_surf_manual$se)
finalTable_avg_surf_manual$low_corr <-  with(finalTable_avg_surf_manual, finalTable_avg_surf_manual$mean - finalTable_avg_surf_manual$se)

ggplot(finalTable_avg_surf_manual, aes(x=time, y=mean, group = treatment, color = treatment)) +
  geom_line() + 
  guides(colour=FALSE)  + 
  geom_ribbon(aes(ymin = low_corr, ymax = high_corr, colour=treatment, group=treatment, fill = treatment ), alpha=.3) +
  # expand_limits(x = 0, y = 0.9) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(limits = c(0.9, 2.6), breaks = seq(0.9, 2.6, by = 0.1)) +
  xlab("Time (s)") + 
  ylab("Norm. fluorescence intensity (A.U.)") + 
  ggtitle("Avg. Surf Norm") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

finalTable_avg_peak_manual <- table.signal_manual %>% group_by(treatment, time) %>% summarize(mean=mean(peak_norm), N = length(peak_norm), sd = sd(peak_norm), se = sd / sqrt(N))

finalTable_avg_peak_manual$high_corr <- with(finalTable_avg_peak_manual, finalTable_avg_peak_manual$mean + finalTable_avg_peak_manual$se)
finalTable_avg_peak_manual$low_corr <-  with(finalTable_avg_peak_manual, finalTable_avg_peak_manual$mean - finalTable_avg_peak_manual$se)

ggplot(finalTable_avg_peak_manual, aes(x=time, y=mean, group = treatment, color = treatment)) +
  geom_line() + 
  guides(colour=FALSE)  + 
  geom_ribbon(aes(ymin = low_corr, ymax = high_corr, colour=treatment, group=treatment, fill = treatment ), alpha=.3) +
  # expand_limits(x = 0, y = 0.9) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  # scale_y_continuous(limits = c(0.9, 1.7), breaks = seq(0.9, 1.7, by = 0.1)) +
  xlab("Time (s)") + 
  ylab("Norm. fluorescence intensity (A.U.)") + 
  ggtitle("Avg. Peak Norm") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

# test fittings
table.signal_manual_treat <- table.signal_manual %>% filter(day == "160525" & treatment == "pN-Blebb")

ggplot(table.signal_manual_treat, aes(x=time, y=peak_norm, group=treatment)) +
  geom_line() + 
  guides(colour=TRUE) + 
  # geom_ribbon(aes(ymin = low_corr, ymax = high_corr, colour=treatment, group=treatment, fill = treatment ), alpha=.3) +
  # expand_limits(x = 0, y = 0.9) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  # scale_y_continuous(limits = c(0.9, 1.7), breaks = seq(0.9, 1.7, by = 0.1)) +
  xlab("Time (s)") + 
  ylab("Norm. fluorescence intensity (A.U.)") + 
  ggtitle("Avg. Peak Norm") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

# reuptake <- subset( table.signal_manual_treat, time >= table.signal_manual_treat$time[ which.max( table.signal_manual_treat$peak_norm ) ] )
reuptake <- subset( table.signal_manual_treat, time >= 14 )

plot(reuptake$time, reuptake$peak_norm)

theta.0 <- min(reuptake$peak_norm) * 0.5  
model.0 <- lm(log(peak_norm - theta.0) ~ time, data=reuptake)  
alpha.0 <- exp(coef(model.0)[1])
beta.0 <- coef(model.0)[2]
start <- list(alpha = alpha.0, beta = beta.0, theta = theta.0)
start

model <- nls(peak_norm ~ alpha * exp(beta * time) + theta , data = reuptake, start = start)

plot(reuptake$time, reuptake$peak_norm )
lines(reuptake$time, predict(model, list(x = reuptake$time)), col = 'skyblue', lwd = 3)

fittingParams <- coef(model)
alpha = fittingParams[1]
beta = fittingParams[2]
theta = fittingParams[3]
y = 1 / exp(1)
y
# y = alpha * exp (beta * time) + theta
# y - theta = alpha * exp(beta * time)
# ( y - theta ) / alpha = exp(beta * time)
# log( (y - theta) / alpha ) = beta * time
time = ( log( (y - theta) / alpha )  / beta )
time

qplot(reuptake$time, reuptake$peak_norm, data =  augment(model)) + 
  geom_hline(yintercept=y) +
  geom_vline(xintercept=time) +
  geom_line(aes(y = .fitted) 
  )

fit <- nls(peak_norm ~ SSasymp(time, yf, y0, log_alpha), data = reuptake)

fittingParams <- coef(fit)
alpha = exp(fittingParams[3])
yf = fittingParams[1]
y0 = fittingParams[2]
y = 1 / exp(1)

tau = log( (y - yf) / (y0 - yf) ) / - alpha
tau

qplot(reuptake$time, reuptake$peak_norm, data =  augment(fit)) + 
  geom_hline(yintercept=y) +
  geom_vline(xintercept=tau) +
  geom_line(aes(y = .fitted) 
  )

# ==============================================================================
# compute peaks
peaks_manual <- table.signal_manual %>% group_by(day, treatment) %>% summarise(max = max(surf_norm))
peaks_manual$deltaMax <- peaks_manual$max - 1

manual_peak_sd <- peaks_manual %>% group_by(treatment) %>% summarize(mean=mean(deltaMax), N = length(deltaMax), sd = sd(deltaMax), se = sd / sqrt(N))
manual_peak_sd

res_peaks_manual <- t.test(deltaMax ~ treatment, data = peaks_manual, paired = TRUE)
res_peaks_manual


# plot peak difference
ggplot(data=peaks_manual, aes(x=treatment, y=deltaMax)) +
  geom_boxplot(outlier.size = 0, outlier.shape = 1) +
  stat_boxplot(geom = 'errorbar', width = 0.2) +
  geom_jitter(width = 0.1) +
  ylab("delta F (exocytosis)") + 
  expand_limits(y = 0) +
  scale_y_continuous(limits = c(0, 1.8), breaks = seq(0, 1.8, by = 0.2)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

# compute taus
table.signal_manual$name <- paste0(table.signal_manual$day, "_", table.signal_manual$treatment)
tau_manual <- calcTau(table.signal_manual, 14)
tau_manual <- tau_manual %>% separate(name, sep ="_", c("day", "treatment"), remove=FALSE)

res_tau_manual <- t.test(tau ~ treatment, data = tau_manual, paired = FALSE)
res_tau_manual 

# box plot of tau values
ggplot(data=tau_manual, aes(x=treatment, y=tau)) +
  geom_boxplot(outlier.size = 0, outlier.shape = 1) +
  stat_boxplot(geom = 'errorbar', width = 0.1) +
  geom_jitter(width = 0.2) +
  ylab("tau") + 
  # expand_limits(y = 20) +
  scale_y_continuous(limits = c(25, 125, breaks = seq(0, 125, by = 10)))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

