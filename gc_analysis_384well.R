# 
# analyzing growth curver data from MIC analysis experiments
# 

{
  library(tidyverse)
  library(ggplot2)
  library(growthcurver)
  library(drc)
  library(dplyr)
  library(scales)
  library(patchwork)
}

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

colors <- c('#b486d2', '#a06dc1', '#844ca9', '#68308D', '#521b76', '#420d64',
            '#fbc068', '#f6ab3b', '#f4a732', '#F09814', '#f89500', '#d98200')

plot_theme <- theme(plot.title = element_text(hjust = 0.5),
                    panel.grid.major=element_line(colour="white"),
                    panel.grid.minor=element_line(colour="white"),
                    legend.key = element_rect(fill = "white"),
                    panel.background = element_rect(fill = "white"),
                    axis.title = element_text(size = 13))

date <- "08062023"
cycle <- "C5"

# before proceeding, annotate the csv file. then, use the following lines to read in the annotated file.
annotated_data <- read.csv(file = paste("../data/mic_", date, "/new_annotated_data.csv", sep = ""))
# annotated_data$Time <- round(as.numeric(seconds(annotated_data$Time)) / 3600, digits = 3) # modifying time to be in hours if needed

# separates out the blank data
blanks <-  startsWith(names(annotated_data), prefix = "blank")
blank_data <- cbind(annotated_data$Time, annotated_data[, blanks, drop = F])

# removes all blank data that has an OD600 above a specified threshold
threshold <- 0.11
blank_data_fixed <- blank_data %>% pivot_longer(!names(blank_data[1]), names_to = "blank", values_to = "OD_600") %>%
  subset(OD_600 < threshold)
names(blank_data_fixed)[1] <- "Time"

# reset the blank data to have problematic measurements removed
blank_data <- blank_data_fixed %>% pivot_wider(names_from = 'blank', values_from = 'OD_600')

# separates sample data from blank data. subtracts average blank OD from sample OD at every time point
sample_data <- cbind(annotated_data[, !blanks, drop = F])

# get row means of all blank data
blank_avg <- rowMeans(blank_data[, c(2:ncol(blank_data))], na.rm = T)

# subtract off the blanks from the sample data
sample_data[, 2:ncol(sample_data)] <- sample_data[, 2:ncol(sample_data)] - blank_avg

######################################################################
# skip these lines if all sample columns are fine
# remove problematic wells from the data set
# problem_wells <- c('fas_L4_C2_1_R2', 'fas_L5_C2_1_R1', 'fas_L5_C2_2.5_R1')
# sample_data <- sample_data[, !names(sample_data) %in% problem_wells]
######################################################################

# calculates growth curve data and splits info about data into individual columns (strain, line, cycle, antibiotic concentration)
{
gc_output <- SummarizeGrowthByPlate(plate = sample_data)
gc_sample_info <- strsplit(gc_output$sample, "_")
gc_output$strain <- unlist(lapply(gc_sample_info, "[", 1))
gc_output$line <- unlist(lapply(gc_sample_info, "[", 2))
gc_output$cycle <- unlist(lapply(gc_sample_info, "[", 3))
gc_output$abx_concentration <- as.numeric(unlist(lapply(gc_sample_info, "[", 4)))
gc_output$replicate <- unlist(lapply(gc_sample_info, "[", 5))
gc_output$label <- as.factor(paste(gc_output$strain, "_", gc_output$line, "_", gc_output$cycle, "_", gc_output$replicate, sep = ""))
}

# creates a list of data frames with relative growth of each isolate at every abx concentration
isolate_labels <- unique(gc_output$label)
isolate_rg <- list()
gc_data <- tibble()
  
for(n in 1:length(isolate_labels)) {
  isolate_gc_data <- gc_output %>% filter(label == isolate_labels[n])
  gc_data <- do.call("rbind", list(gc_data, isolate_gc_data))
}

# summarizes information about area under curve, carrying capacity, and growth rate across all replicates
gc_data_means_auc <- gc_data %>%
  group_by(strain, abx_concentration, line, cycle) %>%
  summarize(average_auc = mean(auc_e), sd_auc = sd(auc_e))

gc_data_means_cc <- gc_data %>%
  group_by(strain, abx_concentration, line, cycle) %>%
  summarize(average_cc = mean(k), sd_cc = sd(k))

gc_data_means_rate <- gc_data %>%
  group_by(strain, abx_concentration, line, cycle) %>%
  summarize(average_rate = mean(r), sd_rate = sd(r))

gc_data_means_auc$label <- as.factor(paste(gc_data_means_auc$strain, "_", gc_data_means_auc$line, sep = ""))
gc_data_means_cc$label <- as.factor(paste(gc_data_means_cc$strain, "_", gc_data_means_cc$line, sep = ""))
gc_data_means_rate$label <- as.factor(paste(gc_data_means_rate$strain, "_", gc_data_means_rate$line, sep = ""))

# creates bar plots with average auc, k, and r grouped by line and antibiotic concentration
auc_plot <- ggplot(data = gc_data_means_auc, aes(x = label, y = average_auc, group = abx_concentration, fill = label)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = average_auc - sd_auc, ymax = average_auc + sd_auc), linewidth = 0.25) +
  plot_theme + 
  scale_fill_manual(values = colors) + 
  ggtitle(paste('Growth Curve Parameters of WT and fas Msm Strains Treated with Isoniazid, ', cycle, sep = '')) + ylab('average area under curve') +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = 'none') +
  facet_wrap(~abx_concentration, nrow = 1)

cc_plot <- ggplot(data = gc_data_means_cc, aes(x = label, y = average_cc, group = abx_concentration, fill = label)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = average_cc - sd_cc, ymax = average_cc + sd_cc), linewidth = 0.25) +
  plot_theme + 
  scale_fill_manual(values = colors) +
  labs(fill = "sample") + ylab('average carrying capacity') +
  ylim(0, 1.5) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.key.height = unit(0.4, 'cm'), 
        legend.title = element_text(size=12), legend.text = element_text(size=12)) +
  facet_wrap(~abx_concentration, nrow = 1)

rate_plot <- ggplot(data = gc_data_means_rate, aes(x = label, y = average_rate, group = abx_concentration, fill = label)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = average_rate - sd_rate, ymax = average_rate + sd_rate), linewidth = 0.25) +
  plot_theme + 
  scale_fill_manual(values = colors) + 
  xlab('isoniazid concentration (µg/mL)') + ylab('average growth rate') +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = 'none') +
  facet_wrap(~abx_concentration, nrow = 1)

pdf(file = paste('../figures/mic-', date,  '/gc_parameters_', cycle, '.pdf', sep = ''), height = 8, width = 12)
wrap_plots(list(auc_plot, cc_plot, rate_plot), nrow = 3)
dev.off()

# saves the data as csv files for each growth parameter
gc_data_means_auc$full.label <- as.factor(paste(gc_data_means_auc$strain, "_", gc_data_means_auc$line, "_", gc_data_means_auc$abx_concentration, "_", gc_data_means_auc$cycle, sep = ""))
gc_data_means_cc$full.label <- as.factor(paste(gc_data_means_cc$strain, "_", gc_data_means_cc$line, "_", gc_data_means_cc$abx_concentration, "_", gc_data_means_cc$cycle, sep = ""))
gc_data_means_rate$full.label <- as.factor(paste(gc_data_means_rate$strain, "_", gc_data_means_rate$line, "_", gc_data_means_rate$abx_concentration, "_", gc_data_means_rate$cycle, sep = ""))

write_csv(gc_data_means_auc, file = paste("../data/mic_", date, "/auc_gc_data-", date, ".csv", sep = ""))
write_csv(gc_data_means_cc, file = paste("../data/mic_", date, "/carrycap_gc_data-", date, ".csv", sep = ""))
write_csv(gc_data_means_rate, file = paste("../data/mic_", date, "/rate_gc_data-", date, ".csv", sep = ""))
