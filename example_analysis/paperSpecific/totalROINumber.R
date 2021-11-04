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
outdir = "/data1/FMP_Docs/Projects/Publication_SynActJ/pHluorinJ_Data/ComparisonROI_Manual_Auto/"

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

sumROI <- manual %>% group_by(treatment) %>% drop_na(NumberROI)  %>% dplyr::summarize(sumManual = sum(NumberROI))

# NOTE need to check if there is a correction for background ROIs
sumROI

# ==============================================================================
# total ROIs after new filter
# ==============================================================================
# further settings
labelSignal = "Spot"
labelBackground = "background"

# Filter settings
sd_multiplicator = 2
peak_filter = 26

table.signal <- read_csv(paste0(indir,"_RawSignal.csv"))
table.background <- read_csv(paste0(indir,"_RawBackground.csv"))

# extracting experimental information from file name
table.signal <- table.signal %>% separate(name, 
                                          sep ="_", c("day", "treatment", "number"), 
                                          remove=FALSE)
table.background <- table.background %>% separate(name, 
                                                  sep ="_", c("day", "treatment", "number"), 
                                                  remove=FALSE)

# filter traces where peak is in the stimulation range (10s - 20s)
table.signal_mean_filter <- subset(table.signal, variable == "mean")

# ------------------------------------------------------------------------------
# computes standard deviation of background overall all traces
table.background_sd <- subset(table.background, variable == "mean")

table.background_sd <- table.background_sd %>% group_by(name, day, treatment, number, roi) %>% dplyr::summarize(sd = sd(value))
table.background_sd$sd_mult <- table.background_sd$sd * sd_multiplicator

# ------------------------------------------------------------------------------
# compute for each trace the average value from the first 4 frames / 2-6 sec
table.signal_mean_before <- table.signal_mean_filter %>% filter(time >= 6 & time <= 8)
table.signal_mean_before_avg <- table.signal_mean_before %>% group_by(name, day, treatment, number, roi) %>% dplyr::summarize(before=mean(value))

# compute for each trace the average value from the next 4 frames / 10-18 sec
table.signal_mean_after <- table.signal_mean_filter %>% filter(time >= 14 & time <= 16)
table.signal_mean_after_avg <- table.signal_mean_after %>% group_by(name, day, treatment, number, roi) %>% dplyr::summarize(after=mean(value))

table.signal_mean_before_after <- left_join(table.signal_mean_before_avg, table.signal_mean_after_avg, by = c("name", "roi", "day", "treatment", "number"))

table.signal_mean_before_after$delta <- table.signal_mean_before_after$after - table.signal_mean_before_after$before

# ------------------------------------------------------------------------------
table.signal_mean_before_after <- left_join(table.signal_mean_before_after, table.background_sd, by = c("name", "roi", "day", "treatment", "number"))

filtered_for_sd <- table.signal_mean_before_after %>% filter(delta > sd_mult)

filtered_sd <- filtered_for_sd %>% select(c(-before, -after, -delta, -sd, -sd_mult) )

filtered_signal_sd <- left_join(filtered_sd, table.signal_mean_filter, by = c("name", "roi", "day", "treatment", "number"))

# ------------------------------------------------------------------------------
# filter traces where peak is in the stimulation range
# ------------------------------------------------------------------------------
# compute peak 
peaks <- filtered_signal_sd %>% group_by(name, roi) %>% dplyr::summarize(value = max(value))
peaks_frame <- left_join(peaks, table.signal_mean_filter, by = c("name", "roi", "value"))

# filter for peak occurance 
filtered_peaks <- peaks_frame %>% filter(time <= peak_filter)
filtered_peaks_2 <- filtered_peaks %>% select(c(-value, -frame, -time, -variable, -value) )

filtered_signal_sd_peak <- left_join(filtered_peaks_2, filtered_signal_sd, by = c("name", "roi", "day", "treatment", "number"))

# compute & plot number of ROIs before filter
unfiltered_singleTP <- subset(table.signal, variable == "mean")
unfiltered_singleTP <- subset(unfiltered_singleTP, time == 0)
roiNumber_noFilter <- unfiltered_singleTP %>% group_by(treatment) %>% dplyr::summarize(sumAuto_nofilter = n())

# ------------------------------------------------------------------------------
# compute & plot number of ROIs after filter
filtered_singleTP <- subset(filtered_signal_sd_peak, time == 0)
roiNumber_newFilter <- filtered_singleTP %>% group_by(treatment) %>% dplyr::summarize(sumAuto_filter = n())

# ------------------------------------------------------------------------------
manulROI_number <- manual %>% group_by(treatment) %>% drop_na(NumberROI)  %>% dplyr::summarize(sumManual = sum(NumberROI))

roiNumber_new_2 <- inner_join(roiNumber_noFilter, manulROI_number, by = "treatment")
roiNumber_new_2 
roiNumber_new_2_longer <- pivot_longer(roiNumber_new_2, !treatment)

ggplot(data=roiNumber_new_2_longer, aes(x=name, y=value)) +
  geom_bar(stat="identity") +
  scale_y_continuous(limits = c(0, 7000), breaks = seq(0, 7000, by = 1000)) +
  ylab("ROI Number") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

# ------------------------------------------------------------------------------
drops <- c("X1","name","frame", "time", "variable", "value")
unfiltered_singleTP_new <- unfiltered_singleTP[, !(names(unfiltered_singleTP) %in% drops)]
unfiltered_singleTP_new$name <- paste0(unfiltered_singleTP_new$day, "_", unfiltered_singleTP_new$treatment, "_", unfiltered_singleTP_new$number, "_", unfiltered_singleTP_new$roi)

drops2 <- c("X1.y","X1.x","name","frame", "time", "variable", "value")
filtered_singleTP_new <- filtered_singleTP[, !(names(filtered_singleTP) %in% drops2)]
filtered_singleTP_new$name <- paste0(filtered_singleTP_new$day, "_", filtered_singleTP_new$treatment, "_", filtered_singleTP_new$number, "_", filtered_singleTP_new$roi)

unfiltered_singleTP_new$kept <- unfiltered_singleTP_new$name %in% filtered_singleTP_new$name

write.csv(unfiltered_singleTP_new, paste0(outdir, "filteredROI.csv"))