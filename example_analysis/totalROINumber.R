setwd("/data1/FMP_Docs/Repositories/plugins_FMP/SynActJ_Shiny/")
library(gridExtra)
library(tidyverse)
library(ggplot2)
source("dataProcessing.R")
source("saveData.R")
source("plotData.R")

# ============================================================================
#
#  DESCRIPTION: Extract and plot number of ROIs Manual vs Automatic analysis
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
#               gridExtra: install.packages("gridExtra")
#               tidyverse: install.packages("tidyverse")
#
#      VERSION: 1.0.0
#      CREATED: 2018-05-24
#     REVISION: 2021-10-21
#
# ============================================================================
# where to get the files
indir = "/data1/FMP_Docs/Projects/Publication_SynActJ/pHluorinJ_Data/RevisedAnalysis/"

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

# ============================================================================
# load data from manual analysis
manual <- read_csv(paste0(indir,"DataTable_ROIs.csv"))

sumROI <- manual %>% group_by(treatment) %>% drop_na(NumberROI)  %>% summarise(sumManual = sum(NumberROI))

# ============================================================================
# load data from automatic analysis
auto <- read_csv(paste0(indir,"_RawSignal.csv"))
auto <- auto %>% separate(name, sep ="_", c("day", "treatment", "number"), remove=FALSE)

# filter traces where peak is in the stimulation range (10s - 20s)
table.signal_mean_filter <- subset(auto, variable == "mean")
peaks <- table.signal_mean_filter %>% group_by(name, roi) %>% summarise(value = max(value))
peaks_frame <- left_join(peaks, table.signal_mean_filter, by = c("name", "roi", "value"))

filtered_peaks <- peaks_frame %>% filter(time < 20)
filtered_peaks_2 <- filtered_peaks %>% select(c(-value, -frame, -time, -variable, -value) )

filtered_signal <- left_join(filtered_peaks_2, table.signal_mean_filter, by = c("name", "roi", "day", "treatment", "number"))

# extracting experimental information from file name
singleData_area <- subset(filtered_signal, time == 0)

# compute & plot number of ROIs
roiNumber <- singleData_area %>% group_by(treatment) %>% summarize(sumAuto = n())

roiNumber <- inner_join(sumROI, roiNumber, by = "treatment")

roiNumber_longer <- pivot_longer(roiNumber, !treatment)
head(roiNumber_longer)

ggplot(data=roiNumber_longer, aes(x=name, y=value)) +
  geom_bar(stat="identity") +
  scale_y_continuous(limits = c(0, 6000), breaks = seq(0, 6000, by = 500)) +
  ylab("ROI Number") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
