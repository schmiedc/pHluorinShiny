setwd("/data1/FMP_Docs/Repositories/plugins_FMP/SynActJ_Shiny/")
source("dataProcessing.R")
source("saveData.R")
source("plotData.R")
library(gridExtra)
library(tidyverse)
library(ggplot2)

# ==============================================================================
#
#  DESCRIPTION: Plot Ctrl vs treatment from output data
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
#               plyr: install.packages("plyr")
#               gridExtra: install.packages("gridExtra")
#               tidyverse: install.packages("tidyverse")
#
#      VERSION: 2.0.0
#      CREATED: 2018-05-24
#     REVISION: 2021-10-21
#
# ==============================================================================
# where to get the files
indir = "/data1/FMP_Docs/Projects/Publication_SynActJ/DataAnalysis/ForZenodo/output/"

# where to save the data
outdir = "/data1/FMP_Docs/Projects/Publication_SynActJ/DataAnalysis/ForZenodo/Routput/"

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

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
# get raw data
table.signal <- collectList(indir, labelSignal, timeResolution)
table.background <- collectList(indir, labelBackground, timeResolution)

# extracting experimental information from file name
table.signal <- table.signal %>% separate(name, 
                                          sep ="_", c("treatment", "number"), 
                                          remove=FALSE)

table.background <- table.background %>% separate(name, 
                                                  sep ="_", c("treatment", "number"), 
                                                  remove=FALSE)

# ------------------------------------------------------------------------------
# Extract and plot number and area of ROIs
# ------------------------------------------------------------------------------
# reduce data for number of ROIs and area
singleData_area <- subset(table.signal, variable == "area")
singleData_area <- subset(singleData_area, time == 0)

# compute & plot number of ROIs
roiNumber <- singleData_area %>% group_by(treatment) %>% dplyr::summarize(count = n())

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
  expand_limits(y = 0) +
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, by = 1)) +
  ylab("Area (Micron)") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

# ------------------------------------------------------------------------------
# Filter extracted traces
# ------------------------------------------------------------------------------
table.signal_mean_filter <- subset(table.signal, variable == "mean")
peaks <- table.signal_mean_filter %>% group_by(name, roi) %>% dplyr::summarise(value = max(value))
peaks_frame <- left_join(peaks, table.signal_mean_filter, by = c("name", "roi", "value"))

# filter traces where peak is in the stimulation range ( < 20s)
filtered_peaks <- peaks_frame %>% filter(time < 20)
filtered_peaks_2 <- filtered_peaks %>% select(c(-value, -frame, -time, -variable, -value) )
filtered_signal <- left_join(filtered_peaks_2, table.signal_mean_filter, by = c("name", "roi", "treatment", "number"))

# ------------------------------------------------------------------------------
# Averaging and normalization
# ------------------------------------------------------------------------------
# calculates average mean intensity per frame per experiment
table.signal_mean <- subset(filtered_signal, variable == "mean")
table.signal_avg <- table.signal_mean %>% group_by(treatment, frame, time) %>% dplyr::summarize(mean=mean(value), N = length(value), sd = sd(value), se = sd / sqrt(N))

table.background_mean <- subset(table.background, variable == "mean")
table.background_avg <- table.background_mean %>% group_by(treatment, frame, time) %>% dplyr::summarize(mean=mean(value), N = length(value), sd = sd(value), se = sd / sqrt(N))

# generate final table
forBackgroundSubtraction <- merge(table.signal_avg, table.background_avg, by=c("treatment", "frame", "time"), suffixes=c(".sig",".back"))

# normalize mean signal with mean background intensity
forBackgroundSubtraction$mean.corr <- forBackgroundSubtraction$mean.sig - forBackgroundSubtraction$mean.back
forBackgroundSubtraction$name <- paste0(forBackgroundSubtraction$treatment)

# surface normalization
surfaceNormalized <- surfaceNormalisation(forBackgroundSubtraction, frameStimulation)

# peak normalization
peakNormalized <- peakNormalisation(surfaceNormalized)

finalTable <- sortFrames(peakNormalized)

# ------------------------------------------------------------------------------
# Plot per individual movie surface and peak normalized data
# ------------------------------------------------------------------------------
ggplot(data=finalTable, aes(x=time, y=surf_norm, group = name, color = name)) +
  geom_line() + 
  theme_light() +
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
  xlab("Time (s)") + 
  ylab("Norm. fluorescence intensity (A.U.)") + 
  ggtitle("Peak Norm") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

# ------------------------------------------------------------------------------
# Plot per treatment surface and peak normalized data
# ------------------------------------------------------------------------------
finalTable_avg_surf <- finalTable %>% group_by(treatment, frame, time) %>% dplyr::summarize(mean=mean(surf_norm), N = length(surf_norm), sd = sd(surf_norm), se = sd / sqrt(N))

ggplot(finalTable_avg_surf, aes(x=time, y=mean, group = treatment, color = treatment)) +
  geom_line() + 
  theme_light() +
  xlab("Time (s)") + 
  ylab("Norm. fluorescence intensity (A.U.)") + 
  ggtitle("Avg. Surf Norm") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))


finalTable_avg_peak <- finalTable %>% group_by(treatment, frame, time) %>% dplyr::summarize(mean=mean(peak_norm), N = length(peak_norm), sd = sd(peak_norm), se = sd / sqrt(N))

ggplot(finalTable_avg_peak, aes(x=time, y=mean, group = treatment, color = treatment)) +
  geom_line() + 
  theme_light() +
  xlab("Time (s)") + 
  ylab("Norm. fluorescence intensity (A.U.)") + 
  ggtitle("Avg. Peak Norm") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

# ------------------------------------------------------------------------------
# Compute and plot peaks based on surface normalization
# ------------------------------------------------------------------------------
# compute peaks
peaks <- finalTable %>% group_by(name) %>% dplyr::summarise(max = max(surf_norm))
peaks$deltaMax <- peaks$max - 1
peaks <- peaks %>% separate(name, sep ="_", c("treatment"), remove=FALSE)
peaks
