setwd("/data1/FMP_Docs/Repositories/plugins_FMP/SynActJ_Shiny/")
source("dataProcessing.R")
source("saveData.R")
source("plotData.R")
library(gridExtra)
library(tidyverse)
library(ggplot2)

# ============================================================================
#
#  DESCRIPTION: Data analysis for a single experiment
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
indir = "/data1/FMP_Docs/Projects/Publication_SynActJ/DataAnalysis/pHluorin_singleMovie/output_160525/"

# where to save the data
outdir = "/data1/FMP_Docs/Projects/Publication_SynActJ/DataAnalysis/pHluorin_singleMovie/Routput/"

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
table.signal <- collectList(indir, labelSignal, timeResolution)
table.background <- collectList(indir, labelBackground, timeResolution)

# ==============================================================================
# create plots for a single data point

# DMSO_10 or DMSO_2
dataName = "DMSO_2"

# compute area for single experiment
singleData_area <- subset(table.signal, variable == "area")
singleData_area <- subset(singleData_area, name == dataName)
singleData_area <- subset(singleData_area, time == 0)

ggplot(data=singleData_area , aes(x=name, y=value)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) +
  expand_limits(y = 0) +
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, by = 1)) +
  ylab("Area (Micron)") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

# compute roi count for single experiment
data <- data.frame(
  name=c("DMSO_2") ,  
  value=c(nrow(singleData_area))
)

ggplot(data=data, aes(x=name, y=value)) +
  geom_bar(stat="identity") +
  ylab("Count") + 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

# raw traces
# plotRawMean_single(table.signal, dataName)
singleData_rawSignal <- subset(table.signal, variable == "mean")
singleData_rawSignal <- subset(singleData_rawSignal, name == dataName)

ggplot(data=singleData_rawSignal, aes(x=time, y=value, group=roi, color = roi)) +
  geom_line() + 
  guides(colour="none")  + 
  theme_light() +
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(limits = c(0, 700), breaks = seq(0, 700, by = 100)) +
  xlab("Time (s)") + 
  ylab("Fluorescence intensity (A.U.)") + 
  ggtitle(paste0("Raw data ", dataName)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

# plotRawMean_single(table.background, dataName)
singleData_rawBack <- subset(table.background, variable == "mean")
singleData_rawBack <- subset(singleData_rawBack, name == dataName)

ggplot(data=singleData_rawBack, aes(x=time, y=value, group=roi, color = roi)) +
  geom_line() + 
  guides(colour="none")  + 
  theme_light() +
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(limits = c(0, 700), breaks = seq(0, 700, by = 100)) +
  xlab("Time (s)") + 
  ylab("Fluorescence intensity (A.U.)") + 
  ggtitle(paste0("Raw data ", dataName)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))


# ==============================================================================
# calculates average mean intensity per frame per experiment
table.signal_mean <- subset(table.signal, variable == "mean")
table.signal_avg <- table.signal_mean %>% group_by(name, frame, time) %>% dplyr::summarize(mean=mean(value), N = length(value), sd = sd(value), se = sd / sqrt(N))

table.background_mean <- subset(table.background, variable == "mean")
table.background_avg <- table.background_mean %>% group_by(name, frame, time) %>% dplyr::summarize(mean=mean(value), N = length(value), sd = sd(value), se = sd / sqrt(N))

# merge averages
forBackgroundSubtraction <- merge(table.signal_avg, table.background_avg, by=c("name", "frame", "time"), suffixes=c(".sig",".back"))

# normalize mean signal with mean background intensity
forBackgroundSubtraction$mean.corr <- forBackgroundSubtraction$mean.sig - forBackgroundSubtraction$mean.back

# ==============================================================================
# normalizations
surfaceNormalized <- surfaceNormalisation(forBackgroundSubtraction, frameStimulation)

peakNormalized <- peakNormalisation(surfaceNormalized)

finalTable <- sortFrames(peakNormalized)

# ==============================================================================
finalTable_avg_sig <- finalTable %>% group_by(frame, time) %>% dplyr::summarize(mean=mean(mean.sig), N = length(mean.sig), sd = sd(mean.sig), se = sd / sqrt(N))

finalTable_avg_sig$high_corr <- with(finalTable_avg_sig, finalTable_avg_sig$mean + finalTable_avg_sig$se)
finalTable_avg_sig$low_corr <- with(finalTable_avg_sig, finalTable_avg_sig$mean - finalTable_avg_sig$se)

ggplot(finalTable_avg_sig, aes(x=time, y=mean)) +
  geom_line() + 
  guides(colour="none")  + 
  geom_ribbon(aes(ymin = low_corr, ymax = high_corr), alpha=.3) +
  # expand_limits(x = 0, y = 0.9) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(limits = c(0, 700), breaks = seq(0, 700, by = 100)) +
  xlab("Time (s)") + 
  ylab("Fluorescence intensity (A.U.)") + 
  ggtitle("Avg. Signal") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

# processed data
finalTable_avg_sig_corr <- finalTable %>% group_by(frame, time) %>% dplyr::summarize(mean=mean(mean.corr), N = length(mean.corr), sd = sd(mean.corr), se = sd / sqrt(N))

finalTable_avg_sig_corr$high_corr <- with(finalTable_avg_sig_corr, finalTable_avg_sig_corr$mean + finalTable_avg_sig_corr$se)
finalTable_avg_sig_corr$low_corr <- with(finalTable_avg_sig_corr, finalTable_avg_sig_corr$mean - finalTable_avg_sig_corr$se)

ggplot(finalTable_avg_sig_corr, aes(x=time, y=mean)) +
  geom_line() + 
  guides(colour="none")  + 
  geom_ribbon(aes(ymin = low_corr, ymax = high_corr), alpha=.3) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(limits = c(0, 700), breaks = seq(0, 700, by = 100)) +
  xlab("Time (s)") + 
  ylab("Corr. Fluorescence intensity (A.U.)") + 
  ggtitle("Avg. Surf Norm") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
# ==============================================================================
# processed data
finalTable_avg_back <- finalTable %>% group_by(frame, time) %>% dplyr::summarize(mean=mean(mean.back), N = length(mean.back), sd = sd(mean.back), se = sd / sqrt(N))

finalTable_avg_back$high_corr <- with(finalTable_avg_back, finalTable_avg_back$mean + finalTable_avg_back$se)
finalTable_avg_back$low_corr <- with(finalTable_avg_back, finalTable_avg_back$mean - finalTable_avg_back$se)

ggplot(finalTable_avg_back, aes(x=time, y=mean)) +
  geom_line() + 
  guides(colour="none")  + 
  geom_ribbon(aes(ymin = low_corr, ymax = high_corr), alpha=.3) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(limits = c(0, 700), breaks = seq(0, 700, by = 100)) +
  xlab("Time (s)") + 
  ylab("Corr. Fluorescence intensity (A.U.)") + 
  ggtitle("Avg. Surf Norm") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

# ==============================================================================
finalTable_avg_surf <- finalTable %>% group_by(frame, time) %>% dplyr::summarize(mean=mean(surf_norm), N = length(surf_norm), sd = sd(surf_norm), se = sd / sqrt(N))

finalTable_avg_surf$high_corr <- with(finalTable_avg_surf, finalTable_avg_surf$mean + finalTable_avg_surf$se)
finalTable_avg_surf$low_corr <-  with(finalTable_avg_surf, finalTable_avg_surf$mean - finalTable_avg_surf$se)

ggplot(finalTable_avg_surf, aes(x=time, y=mean)) +
  geom_line() + 
  guides(colour="none")  + 
  geom_ribbon(aes(ymin = low_corr, ymax = high_corr), alpha=.3) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(limits = c(0.9, 1.7), breaks = seq(0.9, 1.7, by = 0.1)) +
  xlab("Time (s)") + 
  ylab("Norm. fluorescence intensity (A.U.)") + 
  ggtitle("Avg. Surf Norm") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

# compute peak
finalTable_avg_surf['name'] = 'dmso'
peaks <- finalTable_avg_surf %>% group_by(name) %>%  dplyr::summarize(max = max(mean))
peaks

# ------------------------------------------------------------------------------

finalTable_avg_peak <- finalTable %>% group_by(frame, time) %>% dplyr::summarize(mean=mean(peak_norm), N = length(peak_norm), sd = sd(peak_norm), se = sd / sqrt(N))

finalTable_avg_peak$high_corr <- with(finalTable_avg_peak, finalTable_avg_peak$mean + finalTable_avg_peak$se)
finalTable_avg_peak$low_corr <-  with(finalTable_avg_peak, finalTable_avg_peak$mean - finalTable_avg_peak$se)

ggplot(finalTable_avg_peak, aes(x=time, y=mean)) +
  geom_line() + 
  guides(colour="none")  + 
  geom_ribbon(aes(ymin = low_corr, ymax = high_corr), alpha=.3) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  xlab("Time (s)") + 
  ylab("Norm. fluorescence intensity (A.U.)") + 
  ggtitle("Avg. Peak Norm") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

finalTable_avg_peak['name'] = 'dmso'


finalTable_avg_peak <- finalTable_avg_peak %>% 
  rename(
    peak_norm = mean
  )

head(finalTable_avg_peak)

