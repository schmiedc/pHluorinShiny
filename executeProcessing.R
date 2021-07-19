setwd("/data1/FMP_Docs/Repositories/scripts_FMP/pHluorinShiny/")
library(gridExtra)
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

unique(table.signal$name)


# calculates average mean intensity per frame
avg.signal <- calcMean(table.signal)
avg.background <- calcMean(table.background)

# generate final table
finalTable <- processData(indir, frameStimulation, avg.signal, avg.background)

tau <- calcTau(finalTable)

# save files
# writeToCsv(outdir, resultname, table.signal, table.background, finalTable, tau)
#writeToXlsx(outdir, resultname, table.signal, table.background, finalTable, tau)

# ==============================================================================
# Plotting
plot.list <- plotMeans(avg.signal, avg.background, finalTable)

grid.arrange(grobs = plot.list, ncol = 1, nrow = 1, top = "Processing results")

test_plots <- marrangeGrob(plot.list, ncol = 1, nrow = 1, top = "Processing results")

ggsave(plot = test_plots,
       file=paste0(outdir, .Platform$file.sep, "_Result", ".pdf"), 
       width = 297, 
       height = 210, 
       units = "mm") 


plotRawMean_single(table.signal, "DMSO_1")

mean <- subset(table.signal, variable == "mean" )

singleData <- subset(mean, name == "DMSO_1")

single_plot <- ggplot(data=singleData, aes(x=time, y=value, color = roi, group=roi)) +
        geom_line() + 
        guides(colour=FALSE)  + 
        theme_light() +
        xlab("time (s)") + 
        ylab("Fluorescence intensity (a.u.)") + 
        ggtitle(paste0("Raw data ", "DMSO_1")) 

single_plot

# ==============================================================================
# Split between DMSO and pN-Blebb
# Create plots Raw Signal, Raw Background, 
# Mean, Mean-Background, Mean-Background Norm

raw_signal <- plotRawMean(table.signal)
raw_signal_grids <-  marrangeGrob(raw_signal, ncol = 3, nrow = 4, top = "Raw grey values of active boutons")

ggsave(plot = raw_signal_grids,
       file=paste0(outdir, .Platform$file.sep, "_rawAreaBoutons", ".pdf"), 
       width = 297, 
       height = 210, 
       units = "mm") 

grid.arrange(grobs = raw_signal , ncol = 2, top = "Raw grey values of active boutons")
marrangeGrob(raw_signal, ncol = 2, nrow = 2)

raw_background <- plotRawMean(table.background, "Raw grey values of background")

plotRawMean_single(table.signal, "DMSO_1")

raw_signal <- plotRawMean(table.signal)
raw_signal_grids <-  marrangeGrob(raw_signal, ncol = 3, nrow = 4, top = "Raw grey values of active boutons")

