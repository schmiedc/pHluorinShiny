setwd("/data1/FMP_Docs/Repositories/scripts_FMP/pHluorinShiny/")
library(gridExtra)
library(pracma)
source("dataProcessing.R")
source("saveData.R")
source("plotData.R")
source("fitting.R")

# ============================================================================
#
#  DESCRIPTION: Data analysis for Fíji pHlorin workflow
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
indir = "/data1/FMP_Docs/Projects/Publication_SynapseJ/pHluorinJ_Data/TestSet/Output_160525/"

# where to save the data
outdir = "/data1/FMP_Docs/Projects/Publication_SynapseJ/pHluorinJ_Data/TestSet/Output_160525_R/"

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

# save files
# writeToCsv(outdir, resultname, table.signal, table.background, finalTable, tau)
# writeToXlsx(outdir, resultname, table.signal, table.background, finalTable, tau)
# ==============================================================================
# create plots for a single data point

# DMSO_10 or DMSO_2
dataName = "DMSO_2"

# raw traces
# plotRawMean_single(table.signal, dataName)
singleData_rawSignal <- subset(table.signal, variable == "mean")
singleData_rawSignal <- subset(singleData_rawSignal, name == dataName)

ggplot(data=singleData_rawSignal, aes(x=time, y=value, group=roi)) +
  geom_line() + 
  guides(colour=FALSE)  + 
  theme_light() +
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(limits = c(0, 600), breaks = seq(0, 600, by = 100)) +
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

ggplot(data=singleData_rawBack, aes(x=time, y=value, group=roi)) +
  geom_line() + 
  guides(colour=FALSE)  + 
  theme_light() +
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(limits = c(0, 600), breaks = seq(0, 600, by = 100)) +
  xlab("Time (s)") + 
  ylab("Fluorescence intensity (A.U.)") + 
  ggtitle(paste0("Raw data ", dataName)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

# average traces signal
singleData_avgSignal <- subset(avg.signal, name == dataName)

singleData_avgSignal$high <- with(singleData_avgSignal, singleData_avgSignal$mean + singleData_avgSignal$sd)
singleData_avgSignal$low <-  with(singleData_avgSignal, singleData_avgSignal$mean - singleData_avgSignal$sd)

ggplot(data=singleData_avgSignal, aes(x=time, y=mean)) +
  geom_ribbon(aes(ymin = low, ymax = high, colour=name, group=name, fill = name ), alpha=.3) +
  geom_line(colour = "black") + 
  guides(colour=FALSE)  + 
  theme_light() +
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(limits = c(0, 600), breaks = seq(0, 600, by = 100)) +
  xlab("Time (s)") + 
  ylab("Avg. fluorescence intensity (A.U.)") + 
  ggtitle(paste0("Average signal ", dataName)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

# average traces signal background
singleData_avgBackground <- subset(avg.background, name == dataName)
singleData_avgBackground$high <- with(singleData_avgBackground, singleData_avgBackground$mean + singleData_avgBackground$sd)
singleData_avgBackground$low <-  with(singleData_avgBackground, singleData_avgBackground$mean - singleData_avgBackground$sd)

ggplot(data=singleData_avgBackground, aes(x=time, y=mean)) +
  geom_ribbon(aes(ymin = low, ymax = high, colour=name, group=name, fill = name ), alpha=.3) +
  geom_line() + 
  guides(colour=FALSE)  + 
  theme_light() +
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(limits = c(0, 600), breaks = seq(0, 600, by = 100)) +
  xlab("Time (s)") + 
  ylab("Avg. fluorescence intensity (A.U.)") + 
  ggtitle(paste0("Avg background ", dataName)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

# processed data
singleData_finalTable <- subset(finalTable, name == dataName)

ggplot(data=singleData_finalTable, aes(x=time, y=mean.corr)) +
  geom_line() + 
  guides(colour=FALSE)  + 
  theme_light() +
  expand_limits(x = 0, y = 0) +
  # scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  # scale_y_continuous(limits = c(0, 600), breaks = seq(0, 600, by = 100)) +
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
  # scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
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
  # scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  # scale_y_continuous(limits = c(0, 600), breaks = seq(0, 600, by = 100)) +
  xlab("Time (s)") + 
  ylab("Norm. fluorescence intensity (A.U.)") + 
  ggtitle(paste0("Peak Norm ", dataName)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

# get taus
tau
# ==============================================================================
# plots for one dataset treatment vs ctrl
head(finalTable)

# create average for final table per treatment
finalTable_1 <- finalTable %>% separate(name, sep ="_", c("treatment", "number", "correction"))

#
finalTable_2 <- finalTable_1 %>% group_by(treatment, time) %>% summarise(value = mean(surf_norm))

# compute peaks
peaks <- finalTable_2 %>% summarise(max = max(value))

# plots
ggplot(data=finalTable_2, aes(x=time, y=value, group = treatment, color = treatment)) +
  geom_line() + 
  theme_light() +
  expand_limits(x = 0, y = 0.9) +
  # scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(limits = c(0.9, 1.7), breaks = seq(0.9, 1.7, by = 0.1)) +
  xlab("Time (s)") + 
  ylab("Norm. fluorescence intensity (A.U.)") + 
  ggtitle(paste0("Surf Norm ", dataName)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

finalTable_3 <- finalTable_1 %>% group_by(treatment, time) %>% summarise(value = mean(peak_norm))
head(finalTable_3)

ggplot(finalTable_3, aes(x=time, y=value, group = treatment, color = treatment)) +
  geom_line() + 
  theme_light() +
  # expand_limits(x = 0, y = 0.9) +
  # scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  # scale_y_continuous(limits = c(0.9, 1.7), breaks = seq(0.9, 1.7, by = 0.1)) +
  xlab("Time (s)") + 
  ylab("Norm. fluorescence intensity (A.U.)") + 
  ggtitle(paste0("Surf Norm ", dataName)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

# number of ROIs
head(table.signal)
table.signal_roi <- table.signal %>% group_by(name, roi) %>% distinct(roi) 

table.signal_roi <- table.signal_roi %>% separate(name, sep ="_", c("treatment", "number", "correction"))

table.signal_roi %>% group_by(treatment) %>% count()
# table.signal_roi %>% group_by(name) %>% summarise(sum = sum(n))
