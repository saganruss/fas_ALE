# 
# analyzes MIC data from 384-well analysis experiments for each replicate and with IC99 average across replicates
# 

{
  library(tidyverse)
  library(ggplot2)
  library(growthcurver)
  library(drc)
  library(dplyr)
  library(scales)
}

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

colors <- c('#b486d2', '#a06dc1', '#844ca9', '#68308D', '#521b76', '#420d64',
            '#fbc068', '#f6ab3b', '#f4a732', '#F09814', '#f89500', '#d98200')

plot_theme <- theme(plot.title = element_text(hjust = 0.5),
                    panel.grid.major=element_line(colour="lightgray"),
                    panel.grid.minor=element_line(colour="lightgray"),
                    legend.key = element_rect(fill = "white"),
                    panel.background = element_rect(fill = "white"),
                    axis.title = element_text(size = 15))

# function that inputs a drm curve, a y-value, and four parameters. outputs the x-value at inputted y-value.
get_IC_value <- function(drm_fit, y, b, c, d, e) {
  inside_log <- ((d-c)/(y-c))-1
  return(exp(log(inside_log)/b + log(e)))
}


date <- "08062023"
cycle <- "C5"


# prepares data by reading in txt file with raw data and writing to a csv file 
# make sure to delete the first bunch and last bunch of lines before running the next line
raw_data <- read_tsv(file = paste('../data/mic_', date, '/mic-', date, '.txt', sep = ''), locale = locale(encoding = "latin1"))
raw_data <- raw_data %>% dplyr::select(-`T 600`) # remove temperature column
raw_data$Time <- round(as.numeric(seconds(raw_data$Time)) / 3600, digits = 3) # modifying time to be in hours
write_csv(raw_data, file = paste('../data/mic_', date, '/new_annotated_data.csv', sep = ''))


# before proceeding, annotate the csv file. then, use the following lines to read in the annotated file.
annotated_data <- read.csv(file = paste("../data/mic_", date, "/new_annotated_data.csv", sep = ""))
#annotated_data$Time <- round(as.numeric(seconds(annotated_data$Time)) / 3600, digits = 3) # modifying time to be in hours


# separates out the blank data
blanks <-  startsWith(names(annotated_data), prefix = "blank")
blank_data <- cbind(annotated_data$Time, annotated_data[, blanks, drop = F])

# creates a plot of the OD of the blank data below the specified threshold
blank_data_plotting <- blank_data %>% pivot_longer(!names(blank_data[1]), names_to = "blank", values_to = "OD_600")
names(blank_data_plotting)[1] <- "Time"

# plot the blanks before removing problematic wells
ggplot(data = blank_data_plotting, aes(x = Time, y = OD_600, group = blank, color = blank)) +
  geom_point() +
  geom_line()

# removes all blank data that has an OD600 above a specified threshold
threshold <- 0.11
blank_data_fixed <- blank_data %>% pivot_longer(!names(blank_data[1]), names_to = "blank", values_to = "OD_600") %>%
  subset(OD_600 < threshold)
names(blank_data_fixed)[1] <- "Time"

# plot blank data after removing problematic growth
ggplot(data = blank_data_fixed, aes(x = Time, y = OD_600, group = blank, color = blank)) +
  geom_point() +
  geom_line()

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

# splits info about sample data into individual columns (strain, line, cycle, antibiotic concentration) and plots growth curves
{
sample_data_plotting <- sample_data %>% pivot_longer(!names(sample_data[1]), names_to = "sample", values_to = "OD_600")
sample_info <- strsplit(sample_data_plotting$sample, "_")
sample_data_plotting$strain <- as.factor(unlist(lapply(sample_info, "[", 1)))
sample_data_plotting$line <- unlist(lapply(sample_info, "[", 2))
sample_data_plotting$cycle <- unlist(lapply(sample_info, "[", 3))
sample_data_plotting$abx_concentration <- as.factor(as.numeric(unlist(lapply(sample_info, "[", 4))))
sample_data_plotting$replicate <- unlist(lapply(sample_info, "[", 5))
}

sample_data_means <- sample_data_plotting %>%
  group_by(Time, strain, abx_concentration, line) %>%
  summarize(average_OD = mean(OD_600), sd_OD = sd(OD_600))

sample_data_means$label <- as.factor(paste(sample_data_means$strain, "_", sample_data_means$line, "_", sample_data_means$abx_concentration, sep = ""))
sample_data_means$sub.label <- as.factor(paste(sample_data_means$strain, "_", sample_data_means$line, sep = ""))

pdf(file = paste('../figures/mic-', date,  '/growthcurves_', cycle, '.pdf', sep = ''), height = 6, width = 12)
ggplot(data = sample_data_means, aes(x = Time, y = average_OD, group = label, color = abx_concentration)) +
  ggtitle('Growth Curves of WT and fas Lines at Various Concentrations of Isoniazid') + xlab('Time (hours)') +
  geom_point(size = 0.5) +
  geom_line() +
  geom_errorbar(aes(ymin = average_OD - sd_OD, ymax = average_OD + sd_OD), linewidth = 0.25) +
  plot_theme +
  facet_wrap(~sub.label, nrow = 2)
dev.off()

############################################
# spot checking a few samples
sample_data_plotting$label <- paste(sample_data_plotting$strain, '_', sample_data_plotting$line, '_', sample_data_plotting$cycle, '_', sample_data_plotting$replicate, sep = '')
sample_data_plotting$sub.label <- paste(sample_data_plotting$strain, '_', sample_data_plotting$line, '_', sample_data_plotting$cycle, sep = '')
ggplot(data = sample_data_plotting %>% filter(sub.label == 'wt_L1_C5', abx_concentration == 100), aes(x = Time, y = OD_600, group = sample, color = interaction(replicate, abx_concentration))) +
  geom_point(size = 0.5) +
  geom_line() +
  #geom_errorbar(aes(ymin = average_OD - sd_OD, ymax = average_OD + sd_OD), linewidth = 0.25) +
  plot_theme
############################################

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
  
for(n in 1:length(isolate_labels)) {
  isolate_gc_data <- gc_output %>% filter(label == isolate_labels[n])
  rg <- data.frame(rg = (isolate_gc_data$auc_e / isolate_gc_data$auc_e[1]) * 100, 
                   abx_concentration = isolate_gc_data$abx_concentration,
                   isolate = isolate_gc_data$label)
  isolate_rg[[n]] <- rg
}

rg_data <- do.call(rbind, isolate_rg)

{
sample_info <- strsplit(as.character(rg_data$isolate), "_")
rg_data$strain <- as.factor(unlist(lapply(sample_info, "[", 1)))
rg_data$line <- unlist(lapply(sample_info, "[", 2))
rg_data$cycle <- unlist(lapply(sample_info, "[", 3))
rg_data$sub.label <- paste(rg_data$strain, "_", rg_data$line, "_", rg_data$cycle, sep = "")
}
# the following line can be used to remove outlier relative growth data points
# rg_data <- rg_data %>% filter(!(abx_concentration == 1.0 & isolate == "fas_L5_C5"))


# lists to add isolate drm curve information to:
isolate_names <- unique(rg_data$isolate)
isolate_sub.labels <- unique(rg_data$sub.label)
predicted.drc <- list()
abx_prediction_vector <- seq(0, 1000, by = .1)
IC99_values <- list()
rg_at_IC99 <- list()


# for each isolate: fits drm curve and predicts all points on the curve, calculates the IC99 value, 
# and the corresponding relative growth at the IC99 value
for(i in 1:length(isolate_names)) {
  isolate.data <- rg_data %>% filter(isolate == isolate_names[i])
  drm_fit <- drm(rg~abx_concentration, data = isolate.data, fct = LL.4())
  
  predicted.drc[[i]] <- data.frame(replicate.label = isolate_names[i],
                                   sub.label = str_sub(isolate_names[i], 1, -4),
                                   abx_concentration = as.character(abx_prediction_vector),
                                   pred.growth = predict(drm_fit, newdata = data.frame(abx_prediction_vector)))
  
  min_growth_val <- tail(predicted.drc[[i]], n = 1) %>% pull(pred.growth)
  max_growth_val <- head(predicted.drc[[i]], n = 1) %>% pull(pred.growth)
  IC99_values[[i]] <- get_IC_value(drm_fit, ((max_growth_val - min_growth_val) * .01) + min_growth_val,
                                 drm_fit$coefficients[1],
                                 drm_fit$coefficients[2],
                                 drm_fit$coefficients[3],
                                 drm_fit$coefficients[4])
  
  isolate.ic99 <- as.character(round(IC99_values[[i]], digits = 1))
  rg_at_IC99[[i]] <- predicted.drc[[i]] %>% filter(abx_concentration == isolate.ic99) %>% pull(pred.growth)

}

{
IC99_values <- lapply(IC99_values, signif, 3)
predicted_rg <- do.call(rbind, predicted.drc)
ic99s <- data.frame(label = unlist(isolate_names),
                    ic99 = unlist(IC99_values),
                    rg_at_IC99 = unlist(rg_at_IC99))
replicate_info <- strsplit(as.character(ic99s$label), "_")
ic99s$strain <- as.factor(unlist(lapply(replicate_info, "[", 1)))
ic99s$line <- unlist(lapply(replicate_info, "[", 2))
ic99s$cycle <- unlist(lapply(replicate_info, "[", 3))
ic99s$sub.label <- paste(ic99s$strain, "_", ic99s$line, "_", ic99s$cycle, sep = "")
}

IC_data_means <- ic99s %>%
  group_by(sub.label) %>%
  summarize(average_IC99 = mean(ic99), sd_OD = sd(ic99))

IC_data_means$label <- paste(IC_data_means$average_IC99, ' µg/mL', sep = '')



# creates plot for each line that overlays relative growth points and drm curve fit for each replicate. also displays average IC99 value across the replicates.
pdf(file = paste('../figures/mic-', date,  '/rep_avg_rg+ic99_', cycle, '.pdf', sep = ''), height = 6, width = 12)
ggplot(data = rg_data, aes(x = as.numeric(abx_concentration), y = rg, color = sub.label)) +
  geom_point(alpha = 0.5) +
  geom_line(data = predicted_rg, aes(x = as.numeric(abx_concentration), y = pred.growth, group = replicate.label)) +
  geom_text(data = IC_data_means, mapping = aes(x = 1, y = 5, label = label, group = sub.label)) +
  geom_point(data = ic99s, shape = 4, size = 2, stroke = 2, mapping = aes(x = ic99, y = rg_at_IC99, group = label)) +
  scale_x_log10(labels = label_comma(drop0trailing = TRUE)) +
  ggtitle(substitute(paste('Relative Growth and Average IC99 of WT and ', italic('fas Msm'), ' Strains Treated with Isoniazid'))) + xlab('isoniazid (µg/mL)') + ylab('relative growth') +
  plot_theme +
  scale_color_manual(values = colors) +
  theme(legend.position = 'none') +
  facet_wrap(~sub.label, nrow = 2)
dev.off()


# creates data frame containing IC99 value for each isolate and writes to csv file.
IC_data <- data.frame(unlist(isolate_names), unlist(IC99_values))
names(IC_data)[1] <- "Isolate.Name" 
names(IC_data)[2] <- "Average IC99.Value"

write_csv(IC_data, file = paste("../data/mic_", date, "/rep_mic_values-", date, ".csv", sep = ""))
