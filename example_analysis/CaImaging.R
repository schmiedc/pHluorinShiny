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
timeResolution = 0.1

# when stimulation happens
# these frames are used for calcuating the mean intensity
# then this value is used for the surface normalization
frameStimulation = 55

# further settings
labelSignal = "Spot"
labelBackground = "background"

# get raw data
table.signal <- read_csv(paste0(indir,"Ca_RawSignal.csv"))
table.background <- read_csv(paste0(indir,"Ca_RawBackground.csv"))

# calculates average mean intensity per frame per experiment
table.signal_mean <- subset(table.signal, variable == "mean")
table.signal_avg <- table.signal_mean %>% group_by(name, frame, time) %>% summarize(mean=mean(value), N = length(value), sd = sd(value), se = sd / sqrt(N))

table.background_mean <- subset(table.background, variable == "mean")
table.background_avg <- table.background_mean %>% group_by(name, frame, time) %>% summarize(mean=mean(value), N = length(value), sd = sd(value), se = sd / sqrt(N))

# generate final table
forBackgroundSubtraction <- merge(table.signal_avg, table.background_avg, by=c("name", "frame", "time"), suffixes=c(".sig",".back"))

# normalize mean signal with mean background intensity
forBackgroundSubtraction$mean.corr <- forBackgroundSubtraction$mean.sig - forBackgroundSubtraction$mean.back

head(table.signal_mean)

# ==============================================================================
# normalizations
surfaceNormalized <- surfaceNormalisation(forBackgroundSubtraction, frameStimulation)

peakNormalized <- peakNormalisation(surfaceNormalized)

finalTable <- sortFrames(peakNormalized)

# ==============================================================================
finalTable_avg_surf <- finalTable %>% group_by(frame, time) %>% summarize(mean=mean(surf_norm), N = length(surf_norm), sd = sd(surf_norm), se = sd / sqrt(N))

finalTable_avg_surf$high_corr <- with(finalTable_avg_surf, finalTable_avg_surf$mean + finalTable_avg_surf$se)
finalTable_avg_surf$low_corr <-  with(finalTable_avg_surf, finalTable_avg_surf$mean - finalTable_avg_surf$se)

ggplot(finalTable_avg_surf, aes(x=time, y=mean)) +
  geom_line() + 
  guides(colour=FALSE)  + 
  expand_limits(x = 0, y = 0) +
  geom_ribbon(aes(ymin = low_corr, ymax = high_corr), alpha=.3) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20)) +
  scale_y_continuous(limits = c(0, 3.5), breaks = seq(0, 3.5, by =0.5)) +
  xlab("Time (s)") + 
  ylab("Norm. fluorescence intensity (A.U.)") + 
  ggtitle("Surface Normalized") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

# ==============================================================================
# manual analysis
table.signal_manual <- read_csv(paste0(indir,"CaImaging_Manual.csv"))
head(table.signal_manual)

table.signal_manual$Mean <- NULL

table.signal_manual_long <- table.signal_manual %>% gather( "image", "value", -X1)

table.signal_manual_avg <- table.signal_manual_long %>% group_by(X1) %>% summarize(mean=mean(value), N = length(value), sd = sd(value), se = sd / sqrt(N))

table.signal_manual_avg$high_corr <- with(table.signal_manual_avg, table.signal_manual_avg$mean + table.signal_manual_avg$se)
table.signal_manual_avg$low_corr <-  with(table.signal_manual_avg, table.signal_manual_avg$mean - table.signal_manual_avg$se)

head(table.signal_manual_avg)

ggplot(table.signal_manual_avg, aes(x=X1, y=mean)) +
  geom_line() + 
  guides(colour=FALSE)  + 
  expand_limits(x = 0, y = 0) +
  geom_ribbon(aes(ymin = low_corr, ymax = high_corr), alpha=.3) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20)) +
  scale_y_continuous(limits = c(0, 3.5), breaks = seq(0, 3.5, by =0.5)) +
  xlab("Time (s)") + 
  ylab("Norm. fluorescence intensity (A.U.)") + 
  ggtitle("Surface Normalized") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
