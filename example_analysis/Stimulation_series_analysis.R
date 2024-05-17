setwd("/home/schmiedc/FMP_Docs/Repositories/plugins_FMP/SynActJ_Shiny")
library(gridExtra)
library(tidyverse)
library(zoo) # needed for smoothing the traces
source("dataProcessing.R")
source("saveData.R")
source("plotData.R")

# ============================================================================
#
#  DESCRIPTION: Data analysis for stimulation series Agata
#              
#       AUTHOR: Christopher Schmied, 
#      CONTACT: schmied@dzne.de
#     INSITUTE: Leibniz-Forschungsinstitut fuer Molekulare Pharmakologie (FMP)
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
#               zoo: install.packages("zoo")
#
#      VERSION: 1.0.0
#      CREATED: 2024
#     REVISION:
#
# ============================================================================
# where to get the files
indir = "/home/schmiedc/FMP_Docs/Projects/For_Agata/pHluorin/Output/"

# where to save the data
outdir = "/home/schmiedc/FMP_Docs/Projects/For_Agata/pHluorin/TestDev_out/"

# ============================================================================
resultname = "Test"

# Time resolution in seconds
timeResolution = 0.5

# when stimulation happens that was used for
# these frames are used for calcuating the mean intensity
# then this value is used for the surface normalization
frameStimulation = 335

# Setting for smoothing the traces
window_size = 30

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

# get only the mean intensity of the spots
table.signal.mean <- subset(table.signal, variable == "mean" )

# ==============================================================================
# Extract stimulation features
# this list removed stim frame 5, 35 and 455 - most problematic
# stim frame 335 is the big stimulation used for detection of the synapse
stimulation_list_filtered <- list(5,35,65,95,125,155,185,215,245,275,335)

# just get average values for each bin before and after stimulation
# plot them next to each other
# before stimulation is like 3 frames before and including stim frame
# after stim frame
# behavior of a "good" stimulation event:
# stim frame
# stim frame + 1: signal about the same as during stim frame
# stim frame + 2: sharp rise in signal - peak with >10 A.U. of signal
# stim frame + 3 to 18: gradual drop off 
# before 64, 65
# after: 67, 68
# range_before includes the stim_frame
before_stim_range = 3
# range_after is from stim_frame + 2
after_stim_range = 3
# stimulation threshold (A.U.)
stim_threshold = 2

# ------------------------------------------------------------------------------
# computes mean before stimulation
mean_before_stim <- function(data_frame, stim_frame, range){
  
  before_stim_range_lower = stim_frame - ( range - 1 )
  before_stim_range_upper = stim_frame
  return(data_frame %>% filter(frame %in% (before_stim_range_lower:before_stim_range_upper)) %>% summarise( mean = mean(value)))
  
}

# computes mean after stimulation
mean_after_stim <- function(data_frame, stim_frame, range){
  
  after_stim_range_lower = stim_frame + 2
  after_stim_range_upper = after_stim_range_lower + ( range - 1 )
  return(subset_name_roi %>% filter(frame %in% (after_stim_range_lower:after_stim_range_upper)) %>% summarise( mean = mean(value)))
  
}

# Initialize lists to store data
experiment_list <- list()
roi_name_list <- list()
stimulation_name_list <- list()
before_list <- list()
after_list <- list()
difference_list <- list()
responding_stim_list <- list()

# Get unique experiment names
nameTable <- unique(table.signal.mean$name)

list_iterator = 1

for (experimentName in nameTable) {
  
  # Subset data for the current experiment
  subset_name <- subset(table.signal.mean, name == experimentName)
  
  # Get unique ROI names within the current experiment subset
  roiTable <- unique(subset_name$roi)
  
  for (roiName in roiTable) {
    
    # Subset data for the current ROI within the current experiment
    subset_name_roi <- subset(subset_name, roi == roiName)
    
    # Loop over each stimulation frame
    for (stimNumber in stimulation_list_filtered) {
      
      # Compute mean values before and after stimulation
      mean_after_stim_value <- mean_after_stim(subset_name_roi, stimNumber, after_stim_range)
      mean_before_stim_value <- mean_before_stim(subset_name_roi, stimNumber, before_stim_range)
      difference_stimulation <- mean_after_stim_value - mean_before_stim_value 
      
      stim_response = NULL
      
      if (difference_stimulation >= stim_threshold) {
        
        stim_response = TRUE
        
      } else {
        
        stim_response = FALSE
        
      }
      
      experiment_list[[list_iterator]] <- experimentName
      roi_name_list[[list_iterator]] <- roiName
      stimulation_name_list[[list_iterator]] <- stimNumber
      before_list[[list_iterator]] <- mean_before_stim_value
      after_list[[list_iterator]] <- mean_after_stim_value
      difference_list[[list_iterator]] <- difference_stimulation
      responding_stim_list[[list_iterator]] <- stim_response
      
      list_iterator = list_iterator  + 1
      
    }
    
  }
  
}

# Combine lists into a data frame
feature_data <- data.frame(
  name = unlist(experiment_list),
  roi = unlist(roi_name_list),
  stim_frame = unlist(stimulation_name_list),
  before_mean = unlist(before_list),
  after_mean = unlist(after_list),
  difference = unlist(difference_list),
  stim_response = unlist(responding_stim_list)
)

feature_data$stim_frame_char <- as.character(feature_data$stim_frame)

# ==============================================================================
# Filter based on linear thresholding line along start and end points

# user defined settings
# the stim frames used to extract the average after stimulation values
line_start = 65
line_end = 275

# offset to adjust threshold line
offset <- 3
# ------------------------------------------------------------------------------
# Function to calculate y values of a line passing through two points
calculate_slope_intercept <- function(x1, y1, x2, y2) {
  # Calculate the slope (m) of the line
  m <- (y2 - y1) / (x2 - x1)
  
  # Calculate the intercept (b) of the line
  b <- y1 - m * x1
  
  # Return the slope and intercept as a list
  return(list(slope = m, intercept = b))
}

# Function to calculate y values given the slope, intercept, and x values
calculate_y_values <- function(slope, intercept, x_values) {
  # Calculate the corresponding y values using the line equation y = mx + b
  y_values <- slope * x_values + intercept
  
  # Return the y values
  return(y_values)
}

# Initialize lists to store data
filter_experiment_list <- list()
filter_roi_name_list <- list()
filter_slope_positive_list <- list()
filter_stimulation_exceeds_threshold_list <- list()

# Get unique experiment names
nameTable <- unique(table.signal.mean$name)

list_iterator = 1

for (experimentName in nameTable) {
  
  plot.list <- list()
  
  # Subset data for the current experiment
  subset_name <- subset(table.signal.mean, name == experimentName)
  
  # Subset data for the current experiment
  feature_subset_name <- subset(feature_data, name == experimentName)
  
  # Get unique ROI names within the current experiment subset
  roiTable <- unique(subset_name$roi)
  
  for (roiName in roiTable) {
    
    # Subset data for the current ROI within the current experiment
    subset_name_roi <- subset(subset_name, roi == roiName)
    
    # compute a smooth version of the trace
    subset_name_roi$smooth <- rollmean(subset_name_roi$value, k = window_size, fill = NA)

    # Subset data for the current experiment
    feature_subset_name_roi <- subset(feature_subset_name, roi == roiName)
    
    # Define start and end points
    start_value <- feature_subset_name_roi %>% filter(stim_frame == line_start) %>% pull(after_mean)
    end_value <- feature_subset_name_roi %>% filter(stim_frame == line_end) %>% pull(after_mean)
    start_value <- ceiling(start_value + offset)
    end_value <- ceiling(end_value + offset)
    start_value 
    end_value
    
    line_params <- calculate_slope_intercept(line_start, start_value, line_end, end_value)
    
    x_values <- seq(0, nrow(subset_name_roi), length.out = nrow(subset_name_roi))  # Generate x values from 0 to x2 for plotting
    
    # Calculate predicted y values using the returned predicted_y function
    subset_name_roi$threshold_line <- calculate_y_values(line_params$slope, line_params$intercept, x_values)
    
    # get only the window for determining traces to trash
    subset_name_roi_window <- subset_name_roi %>% filter(frame >= line_start, frame <= line_end)
    
    # if smooth line goes over threshold line the trace will be discarded
    slope_positive = FALSE
    stimulation_exceeds_threshold = FALSE
    
    if  ( line_params$slope > 0) {
      
      slope_positive = TRUE
      
    } else {
      
      slope_positive = FALSE
      
      for (frame_number in subset_name_roi_window$frame) {
        
        subset_name_roi_window_frame <- subset(subset_name_roi_window, frame == frame_number)
        
        if  (subset_name_roi_window_frame$smooth > subset_name_roi_window_frame$threshold_line) {
          
          stimulation_exceeds_threshold = TRUE
          break
          
        } else {
          
          stimulation_exceeds_threshold = FALSE
          
        }
        
      }
      
    }
    
    if (slope_positive | stimulation_exceeds_threshold) {
      theme_background = "red"
      
    } else {
      
      theme_background = "lightgreen"
      
    }
    
    plot.list[[roiName]] <- ggplot() +
      geom_line(data = subset_name_roi, aes(x = frame, y = value, color = "cyan"), na.rm=TRUE, show.legend = FALSE) +
      geom_line(data = subset_name_roi, aes(x = frame, y = smooth, color = "magenta"), na.rm=TRUE, show.legend = FALSE) +
      geom_line(data = subset_name_roi, aes(x = frame, y = threshold_line, color = "green"), na.rm=TRUE, show.legend = FALSE ) +
      theme_bw(base_size = 5) +
      theme(plot.background = element_rect(fill = theme_background ) ) + 
      xlab("Frame number") +
      ylab("Fluorescence (A.U.)") +
      ggtitle(paste0(roiName,
                     "\n",
                     "Positive slope: ",
                     slope_positive ,
                     "\n",
                     "Threshold exceeded: ",
                     stimulation_exceeds_threshold ) ) 
    
    filter_experiment_list[[list_iterator]] <- experimentName
    filter_roi_name_list[[list_iterator]] <- roiName
    filter_slope_positive_list[[list_iterator]] <- slope_positive
    filter_stimulation_exceeds_threshold_list[[list_iterator]] <- stimulation_exceeds_threshold
    
    list_iterator = list_iterator + 1
    
  }
  
  test_plots <- marrangeGrob(plot.list, 
                             ncol = 2, 
                             nrow = 3, 
                             top = "Processing results",
                             layout_matrix = matrix(1:12, 4, 3, TRUE) )
  
  
  ggsave(plot = test_plots,
         file=paste0(outdir, .Platform$file.sep, 'Filter_', experimentName, ".pdf"), 
         width = 297, 
         height = 210, 
         units = "mm") 
  
  plot.list <- NULL
  
}

# Combine lists into a data frame
first_feature_filter <- data.frame(
  name = unlist(filter_experiment_list),
  roi = unlist(filter_roi_name_list),
  slope_positive = unlist(filter_slope_positive_list),
  stimulation_exceeds_threshold = unlist(filter_stimulation_exceeds_threshold_list)
)

sum(first_feature_filter$stimulation_exceeds_threshold, na.rm = TRUE)

# ==============================================================================
head(table.signal.mean) # original trace
head(feature_data) # info about features plus responding stimulations
head(first_feature_filter) # information about filter

# filter traces and features
first_feature_filter_removed <- filter(first_feature_filter, !slope_positive)
first_feature_filter_removed2 <- filter(first_feature_filter_removed , !stimulation_exceeds_threshold)

table.signal.mean_filtered <- table.signal.mean %>% semi_join(first_feature_filter_removed2, by = c("name", "roi"))
feature_data_filtered <- feature_data %>% semi_join(first_feature_filter_removed2, by = c("name", "roi"))

nameTable = unique(table.signal.mean_filtered$name)

for (nameName in nameTable) {
  
  plot.list <- list()
  
  name_subset <- subset(table.signal.mean_filtered, name == nameName )
  feature_data_filtered_subset <- subset(feature_data_filtered, name == nameName)
  
  roiTable = unique(name_subset$roi)
  
  for (roiNumber in roiTable) {
    
    roi_subset <- subset(name_subset, roi == roiNumber)
    feature_data_filtered_subset_roi <- subset(feature_data_filtered_subset, roi == roiNumber)

    plot.list[[roiNumber]]  <- ggplot() +
      
      geom_line(data = roi_subset, aes(x=frame, y=value), size = 0.2) +
      
      xlab("Frame") + 
      ylab("Fluorescence (a.u.)") +
      
      # 
      geom_rect(
        aes(xmin = stim_frame, xmax = stim_frame + 30, fill = stim_response), 
        ymin = -Inf, 
        ymax = Inf, 
        alpha = 0.2, 
        data = feature_data_filtered_subset_roi,
        show.legend = FALSE
      ) + 
      
      geom_vline(
        aes(xintercept = as.numeric(stim_frame_char), color = stim_response), 
        data = feature_data_filtered_subset_roi,
        #colour = "blue", 
        linetype = "dotted",
        size = 0.2
      ) + 
      
      theme_bw(base_size = 5) +
      ggtitle(roiNumber)
    
  }

  test_plots <- marrangeGrob(plot.list, 
                             ncol = 3, 
                             nrow = 4, 
                             top = "Processing results",
                             layout_matrix = matrix(1:12, 4, 3, TRUE) )
  
  
  ggsave(plot = test_plots,
         file=paste0(outdir, .Platform$file.sep, "Stim_", nameName, ".pdf"), 
         width = 297, 
         height = 210, 
         units = "mm") 
  
  plot.list <- NULL
  
}

# ==============================================================================
# Plot extracted features

for (nameName in nameTable) {
  
  plot_list_mean_before <- list()
  plot_list_difference <- list()
  
  feature_subset <- subset(feature_data, name == nameName )
  
  roiTable = unique(feature_subset$roi)
  
  for (roiNumber in roiTable) {
    
    feature_subset_roi <- subset(feature_subset, roi == roiNumber )
    
    plot_list_mean_before[[roiNumber]] <- ggplot(data=feature_subset_roi, aes(x=fct_inorder(stim_frame_char), y=before_mean)) +
      geom_bar(stat="identity") +
      xlab("Stimulations frames") + 
      ylab("Fluorescence (a.u.)") +
      ggtitle(roiNumber)
    
    plot_list_difference[[roiNumber]] <- ggplot(data=feature_subset_roi, aes(x=fct_inorder(stim_frame_char), y=difference)) +
      geom_bar(stat="identity") +
      xlab("Stimulations frames") + 
      ylab("After - Before (a.u.)") +
      ggtitle(roiNumber)
    
  }
  
  plot_list_mean_before_grid <- marrangeGrob(plot_list_mean_before, 
                                             ncol = 3, 
                                             nrow = 4, 
                                             top = "Mean signal before stimulation",
                                             layout_matrix = matrix(1:12, 4, 3, TRUE) )
  
  
  ggsave(plot = plot_list_mean_before_grid,
         file=paste0(outdir, .Platform$file.sep, "Mean_before_", nameName, ".pdf"), 
         width = 297, 
         height = 210, 
         units = "mm") 
  
  plot_list_difference_grid <- marrangeGrob(plot_list_difference, 
                                            ncol = 3, 
                                            nrow = 4, 
                                            top = "Mean signal before stimulation",
                                            layout_matrix = matrix(1:12, 4, 3, TRUE) )
  
  
  ggsave(plot = plot_list_difference_grid,
         file=paste0(outdir, .Platform$file.sep, "Difference_", nameName, ".pdf"), 
         width = 297, 
         height = 210, 
         units = "mm") 
  
  plot.list <- NULL
  
}


# ==============================================================================
# Filtering of traces
lots_of_trash = "2312127_1_iNWT0311_VGLUTpHlenti1uL_100nMJFX650_1.3mMCa_10x4AP_t5x30f_200AP_t335_500msframe_2_MMImages.ome"
good_traces = "2312127_1_iNWT0311_VGLUTpHlenti1uL_HaloFBPWT1uL_100nMJFX650_1.3mMCa_10x4AP_t5x30f_200AP_t335_500msframe_3_MMImages.ome"

signal_subset <- subset(table.signal.mean, name == good_traces)
signal_subset_roi <- subset(signal_subset, roi == '010')

feature_subset <- subset(feature_data, name == good_traces)
feature_subset_roi <- subset(feature_subset, roi == '010')

