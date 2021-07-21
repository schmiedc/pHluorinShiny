setwd("/data1/FMP_Docs/Repositories/plugins_FMP/pHluorinShiny/")
library(gridExtra)
source("dataProcessing.R")
source("saveData.R")
source("plotData.R")
source("fitting.R")

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
# indir = "/data1/FMP_Docs/Projects/Publication_SynapseJ/pHluorinJ_Data/AutomaticAnalysisOut/"
indir = "/data1/FMP_Docs/Projects/Publication_SynapseJ/pHluorinJ_Data/TestSet/Output_160525/"

# where to save the data
outdir = indir
# outdir = "/home/schmiedc/Desktop/Output_160525/"

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

# calculates average mean intensity per frame
avg.signal <- calcMean(table.signal)
avg.background <- calcMean(table.background)

# generate final table
finalTable <- processData(indir, frameStimulation, avg.signal, avg.background)

tau <- calcTau(finalTable)
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
  guides(colour=FALSE)  + 
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
  guides(colour=FALSE)  + 
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

# average traces signal
singleData_avgSignal <- subset(avg.signal, name == dataName)

singleData_avgSignal$high <- with(singleData_avgSignal, singleData_avgSignal$mean + singleData_avgSignal$se)
singleData_avgSignal$low <-  with(singleData_avgSignal, singleData_avgSignal$mean - singleData_avgSignal$se)

ggplot(data=singleData_avgSignal, aes(x=time, y=mean)) +
  geom_ribbon(aes(ymin = low, ymax = high, colour=name, group=name, fill = name ), alpha=.3) +
  geom_line(colour = "black") + 
  guides(colour=FALSE)  + 
  theme_light() +
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(limits = c(0, 700), breaks = seq(0, 700, by = 100)) +
  xlab("Time (s)") + 
  ylab("Avg. fluorescence intensity (A.U.)") + 
  ggtitle(paste0("Average signal ", dataName)) +
  theme(panel.grid.major = elemeunt_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

# average traces signal background
singleData_avgBackground <- subset(avg.background, name == dataName)
singleData_avgBackground$high <- with(singleData_avgBackground, singleData_avgBackground$mean + singleData_avgBackground$se)
singleData_avgBackground$low <-  with(singleData_avgBackground, singleData_avgBackground$mean - singleData_avgBackground$se)

ggplot(data=singleData_avgBackground, aes(x=time, y=mean)) +
  geom_ribbon(aes(ymin = low, ymax = high, colour=name, group=name, fill = name ), alpha=.3) +
  geom_line() + 
  guides(colour=FALSE)  + 
  theme_light() +
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(limits = c(0, 700), breaks = seq(0, 700, by = 100)) +
  xlab("Time (s)") + 
  ylab("Avg. fluorescence intensity (A.U.)") + 
  ggtitle(paste0("Avg background ", dataName)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
u
# processed data
singleData_finalTable <- subset(finalTable, name == dataName)
singleData_finalTable$high_corr <- with(singleData_finalTable, singleData_finalTable$mean.corr + singleData_finalTable$sd.sig)
singleData_finalTable$low_corr <-  with(singleData_finalTable, singleData_finalTable$mean.corr - singleData_finalTable$sd.sig)

ggplot(data=singleData_finalTable, aes(x=time, y=mean.corr)) +
  geom_ribbon(aes(ymin = low_corr, ymax = high_corr, colour=name, group=name, fill = name ), alpha=.3) +
  geom_line() + 
  guides(colour=FALSE)  + 
  theme_light() +
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(limits = c(0, 700), breaks = seq(0, 700, by = 100)) +
  xlab("Time (s)") + 
  ylab("Norm. fluorescence intensity (A.U.)") + 
  ggtitle(paste0("Background sub ", dataName)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

ggplot(data=singleData_finalTable, aes(x=time, y=surf_norm)) +
  geom_line() + 
  guides(colour=FALSE)  + 
  theme_light() +
  expand_limits(x = 0, y = 0.9) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(limits = c(0.9, 1.7), breaks = seq(0.9, 1.7, by = 0.1)) +
  xlab("Time (s)") + 
  ylab("Norm. fluorescence intensity (A.U.)") + 
  ggtitle(paste0("Surf Norm ", dataName)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

ggplot(data=singleData_finalTable, aes(x=time, y=peak_norm)) +
  geom_line() + 
  guides(colour=FALSE)  + 
  theme_light() +
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  # scale_y_continuous(limits = c(0, 600), breaks = seq(0, 600, by = 100)) +
  xlab("Time (s)") + 
  ylab("Norm. fluorescence intensity (A.U.)") + 
  ggtitle(paste0("Peak Norm ", dataName)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

# compute peak
head(singleData_finalTable)
peaks <- singleData_finalTable %>% summarise(max = max(surf_norm))
peaks
tau
