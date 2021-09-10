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

# ------------------------------------------------------------------------------
outputDirectory = "/data1/FMP_Docs/Projects/Publication_SynapseJ/pHluorinJ_Data/checkFilter/"
sd_multiplicator = 2
peak_filter = 26

# ==============================================================================
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

# ==============================================================================
# computes standard deviation of background overall all traces
table.background_sd <- subset(table.background, variable == "mean")

table.background_sd <- table.background_sd %>% group_by(name, day, treatment, number, roi) %>% summarize(sd = sd(value))
table.background_sd$sd_mult <- table.background_sd$sd * sd_multiplicator

# ==============================================================================
table.signal_mean_filter <- subset(table.signal, variable == "mean")


raw_signal_unfiltered <- plotRawMean(table.signal_mean_filter)
raw_signal_grids_unfiltered <-  marrangeGrob(raw_signal_unfiltered, ncol = 3, nrow = 4, top = "Raw grey values of filtered traces")
ggsave(plot = raw_signal_grids_unfiltered,
       file=file.path(outputDirectory, paste0(resultname, "_unfiltered", sd_multiplicator, ".pdf") ), 
       width = 297, 
       height = 210, 
       units = "mm") 

# ==============================================================================
# compute for each trace the average value from the first 4 frames / 2-6 sec
table.signal_mean_before <- table.signal_mean_filter %>% filter(time >= 6 & time <= 8)
table.signal_mean_before_avg <- table.signal_mean_before %>% group_by(name, day, treatment, number, roi) %>% summarize(before=mean(value))

# compute for each trace the average value from the next 4 frames / 10-18 sec
table.signal_mean_after <- table.signal_mean_filter %>% filter(time >= 14 & time <= 16)
table.signal_mean_after_avg <- table.signal_mean_after %>% group_by(name, day, treatment, number, roi) %>% summarize(after=mean(value))

table.signal_mean_before_after <- left_join(table.signal_mean_before_avg, table.signal_mean_after_avg, by = c("name", "roi", "day", "treatment", "number"))

table.signal_mean_before_after$delta <- table.signal_mean_before_after$after - table.signal_mean_before_after$before

# ==============================================================================
table.signal_mean_before_after <- left_join(table.signal_mean_before_after, table.background_sd, by = c("name", "roi", "day", "treatment", "number"))

filtered_for_sd <- table.signal_mean_before_after %>% filter(delta > sd_mult)

filtered_sd <- filtered_for_sd %>% select(c(-before, -after, -delta, -sd, -sd_mult) )

filtered_signal_sd <- left_join(filtered_sd, table.signal_mean_filter, by = c("name", "roi", "day", "treatment", "number"))

raw_signal_sd <- plotRawMean(filtered_signal_sd)
raw_signal_grids_sd <-  marrangeGrob(raw_signal_sd, ncol = 3, nrow = 4, top = "Raw grey values of filtered traces")
ggsave(plot = raw_signal_grids_sd,
       file=file.path(outputDirectory, paste0(resultname, "_keep_sd", sd_multiplicator, ".pdf") ), 
       width = 297, 
       height = 210, 
       units = "mm") 
# ------------------------------------------------------------------------------
# save rejected
filtered_for_sd_reject <- table.signal_mean_before_after %>% filter(delta < sd_mult)

filtered_sd_reject <- filtered_for_sd_reject %>% select(c(-before, -after, -delta, -sd, -sd_mult) )

filtered_signal_sd_reject <- left_join(filtered_sd_reject, table.signal_mean_filter, by = c("name", "roi", "day", "treatment", "number"))

raw_signal_sd_reject <- plotRawMean(filtered_signal_sd_reject)
raw_signal_grids_sd_reject <-  marrangeGrob(raw_signal_sd_reject, ncol = 3, nrow = 4, top = "Raw grey values of filtered traces")
ggsave(plot = raw_signal_grids_sd_reject,
       file=file.path(outputDirectory, paste0(resultname, "_reject_sd", sd_multiplicator, ".pdf") ), 
       width = 297, 
       height = 210, 
       units = "mm") 

# ==============================================================================
# filter traces where peak is in the stimulation range (20s)

# compute peak 
peaks <- filtered_signal_sd %>% group_by(name, roi) %>% summarise(value = max(value))
peaks_frame <- left_join(peaks, table.signal_mean_filter, by = c("name", "roi", "value"))

# filter for peak occurance 
filtered_peaks <- peaks_frame %>% filter(time <= peak_filter)
filtered_peaks_2 <- filtered_peaks %>% select(c(-value, -frame, -time, -variable, -value) )

filtered_signal_sd_peak <- left_join(filtered_peaks_2, filtered_signal_sd, by = c("name", "roi", "day", "treatment", "number"))

# create plot grid and save
raw_signal <- plotRawMean(filtered_signal_sd_peak)
raw_signal_grids <-  marrangeGrob(raw_signal, ncol = 3, nrow = 4, top = "Raw grey values of filtered traces")
ggsave(plot = raw_signal_grids,
       file=file.path(outputDirectory, paste0(resultname, "_keep_sd", sd_multiplicator,"_time", peak_filter,".pdf") ), 
       width = 297, 
       height = 210, 
       units = "mm") 

# ------------------------------------------------------------------------------
# save rejected
# filter for peak occurance 
filtered_peaks_reject <- peaks_frame %>% filter(time >= peak_filter)
filtered_peaks_2_reject <- filtered_peaks_reject %>% select(c(-value, -frame, -time, -variable, -value) )

filtered_signal_sd_peak_reject <- left_join(filtered_peaks_2_reject, filtered_signal_sd, by = c("name", "roi", "day", "treatment", "number"))

# create plot grid and save
raw_signal_reject <- plotRawMean(filtered_signal_sd_peak_reject)
raw_signal_grids_reject <-  marrangeGrob(raw_signal_reject, ncol = 3, nrow = 4, top = "Raw grey values of filtered traces")
ggsave(plot = raw_signal_grids_reject,
       file=file.path(outputDirectory, paste0(resultname, "_reject_time", peak_filter,".pdf") ), 
       width = 297, 
       height = 210, 
       units = "mm") 

# ==============================================================================
# calculates average mean intensity per frame per experiment
table.signal_mean <- subset(table.signal, variable == "mean")
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
save_table <- finalTable %>% select(.id, time, peak_norm)

save_table_wider <- save_table %>% spread(.id, peak_norm)

library("xlsx")
write.xlsx(save_table_wider, file.path(outputDirectory, paste0(resultname, "_time30.xlsx")), sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)

# write.xlsx(save_table_wider, file.path(outputDirectory, paste0(resultname, "_sd", sd_multiplicator,"_time", peak_filter, ".xlsx")), sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)
# ==============================================================================
finalTable_avg_surf <- finalTable %>% group_by(treatment, frame, time) %>% summarize(mean=mean(surf_norm), N = length(surf_norm), sd = sd(surf_norm), se = sd / sqrt(N))

finalTable_avg_surf$high_corr <- with(finalTable_avg_surf, finalTable_avg_surf$mean + finalTable_avg_surf$se)
finalTable_avg_surf$low_corr <-  with(finalTable_avg_surf, finalTable_avg_surf$mean - finalTable_avg_surf$se)

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
res_peaks

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
