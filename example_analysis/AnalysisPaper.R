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
indir = "/data1/FMP_Docs/Projects/Publication_SynapseJ/pHluorinJ_Data/AutomaticAnalysisOut/"

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

# ==============================================================================
# plots for one dataset treatment vs ctrl
head(finalTable)

# create average for final table per treatment
finalTable_1 <- finalTable %>% separate(name, sep ="_", c("day", "treatment", "number"), remove=FALSE)
# finalTable_1 <- finalTable %>% separate(name, sep ="_", c("treatment", "number"), remove=FALSE)
head(finalTable_1)
#
finalTable_2 <- finalTable_1 %>% group_by(treatment, time) %>% summarise(value = mean(surf_norm))
head(finalTable_1)
# compute peaks
peaks <- finalTable_2 %>% summarise(max = max(value))

# plots
ggplot(data=finalTable_2, aes(x=time, y=value, group = treatment, color = treatment)) +
  geom_line() + 
  theme_light() +
  expand_limits(x = 0, y = 0.9) +
  # scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(limits = c(0.9, 2.5), breaks = seq(0.9, 2.5, by = 0.1)) +
  xlab("Time (s)") + 
  ylab("Norm. fluorescence intensity (A.U.)") + 
  ggtitle("Surf Norm ") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

finalTable_3 <- subset(finalTable_1, day %in% c("160629", "160525" ) ) %>% group_by(treatment, time) %>% summarise(value = mean(peak_norm))

finalTable_3 <- subset(finalTable_1, day %in% c( "160629" ) ) %>% group_by(treatment, time) %>% summarise(value = mean(surf_norm))
head(finalTable_3)

ggplot(finalTable_3, aes(x=time, y=value, group = treatment, color = treatment)) +
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

head(finalTable_1)
finalTable_4 <- finalTable_1 %>% group_by(day, treatment, time) %>% summarise(value = mean(peak_norm))
head(finalTable_4)

ggplot(subset(finalTable_4, treatment %in% "pN-Blebb"), aes(x=time, y=value, group = day, color = day)) +
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

# number of ROIs
table.signal_roi <- table.signal %>% group_by(name, roi) %>% distinct(roi) 

table.signal_roi <- table.signal_roi %>% separate(name, sep ="_", c("day", "treatment", "number"))

table.signal_roi %>% group_by(treatment) %>% count()
# table.signal_roi %>% group_by(name) %>% summarise(sum = sum(n))
