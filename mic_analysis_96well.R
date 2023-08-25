# analyzing MIC data from 96-well analysis experiments

library(tidyverse)
library(ggplot2)
library(growthcurver)
library(drc)

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

date <- "07212023"
cycle <- "C5"


# prepares data by reading in txt file with annotated data and writes to a csv file 
annotated_data <- read_delim(file = paste("../data/mic_", date, "/annotated_mic-", date, ".txt", sep = ""), delim = '\t')
annotated_data$Time <- round(as.numeric(seconds(annotated_data$Time)) / 3600, digits = 3) # modifying time to be in hours
annotated_data <- annotated_data %>% dplyr::select(-`T 600`) # remove temperature column
write_csv(annotated_data, file = paste("../data/mic_", date, "/annotated_mic-", date, ".csv", sep = ""))


# separates out the blank data
blanks <-  startsWith(names(annotated_data), prefix = "blank")
blank_data <- cbind(annotated_data$Time, annotated_data[, blanks, drop = F])

# removes all blanks that have an OD600 above a specified threshold
threshold <- 0.11
outlier_blank_idx <- c()

for (i in 2:length(blank_data)) {
  if(any(blank_data[i] > threshold)) {
    outlier_blank_idx <- append(outlier_blank_idx, c(i))
  }
}

blank_data <- blank_data[-outlier_blank_idx]


# creates a plot of the OD of the blank data below the specified threshold
blank_data_plotting <- blank_data %>% pivot_longer(!names(blank_data[1]), names_to = "blank", values_to = "OD_600")
names(blank_data_plotting)[1] <- "Time"

ggplot(data = blank_data_plotting, aes(x = Time, y = OD_600, group = blank, color = blank)) +
  geom_point() +
  geom_line()


# separates sample data from blank data. subtracts average blank OD from sample OD at every time point
sample_data <- cbind(annotated_data[, !blanks, drop = F])
sample_data[, 2:ncol(sample_data)] <- sample_data[, 2:ncol(sample_data)] - rowMeans(blank_data[, 2:ncol(blank_data)])


# splits info about sample data into individual columns (strain, line, cycle, antibiotic concentration) and plots growth curves
{
sample_data_plotting <- sample_data %>% pivot_longer(!names(sample_data[1]), names_to = "sample", values_to = "OD_600")
sample_info <- strsplit(sample_data_plotting$sample, "_")
sample_data_plotting$strain <- as.factor(unlist(lapply(sample_info, "[", 1)))
sample_data_plotting$line <- unlist(lapply(sample_info, "[", 2))
sample_data_plotting$cycle <- unlist(lapply(sample_info, "[", 3))
sample_data_plotting$abx_concentration <- unlist(lapply(sample_info, "[", 4))
}

ggplot(data = sample_data_plotting, aes(x = Time, y = OD_600, group = sample, color = sample)) +
  ggtitle('Growth Curves of WT and fas Lines at Various Concentrations of Isoniazid') + xlab('Time (hours)') +
  geom_point() +
  geom_line()

# calculates growth curve data and splits info about data into individual columns (strain, line, cycle, antibiotic concentration)
{
gc_output <- SummarizeGrowthByPlate(plate = sample_data)
gc_sample_info <- strsplit(gc_output$sample, "_")
gc_output$strain <- unlist(lapply(gc_sample_info, "[", 1))
gc_output$line <- unlist(lapply(gc_sample_info, "[", 2))
gc_output$cycle <- unlist(lapply(gc_sample_info, "[", 3))
gc_output$abx_concentration <- unlist(lapply(gc_sample_info, "[", 4))
gc_output$label <- as.factor(paste(gc_output$strain, "_", gc_output$line, "_", gc_output$cycle, sep = ""))
}

# creates a list of data frames with relative growth of each isolate at every abx concentration
isolate_labels <- unique(gc_output$label)
isolate_rg <- list()
  
for(n in 1:length(isolate_labels)) {
  i <- 4
  isolate_gc_data <- gc_output %>% filter(label == isolate_labels[n])
  rg <- data.frame(rg = (isolate_gc_data$auc_e / isolate_gc_data$auc_e[1]) * 100, 
                   abx_concentration = as.numeric(isolate_gc_data$abx_concentration),
                   isolate = isolate_gc_data$label)
  isolate_rg[[n]] <- rg
}

rg_data <- do.call(rbind, isolate_rg)

# the following line can be used to remove outlier relative growth data points
# was used for experiment 1, cycle 5, fas L5 at abx concentration of 1.0:
# rg_data <- rg_data %>% filter(!(abx_concentration == 1.0 & isolate == "fas_L5_C5"))


# lists to add isolate drm curve information to:
isolate_names <- unique(rg_data$isolate)
predicted.drc <- list()
abx_prediction_vector <- seq(0, 1000, by = .1)
IC99_values <- list()
rg_at_IC99 <- list()


# for each isolate: fits drm curve and predicts all points on the curve, calculates the IC99 value, 
# and the corresponding relative growth at the IC99 value
for(i in 1:length(isolate_names)) {
  isolate.data <- rg_data %>% filter(isolate == isolate_names[i])
  drm_fit <- drm(rg~abx_concentration, data = isolate.data, fct = LL.4(fixed = c(NA, NA, NA, NA)))
  
  predicted.drc[[i]] <- data.frame(isolate = isolate_names[i],
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

IC99_values <- lapply(IC99_values, signif, 3)
predicted_rg <- do.call(rbind, predicted.drc)
ic99s <- data.frame(isolate = unlist(isolate_names),
                    ic99 = unlist(IC99_values),
                    rg_at_IC99 = unlist(rg_at_IC99))
ic99s$ic99.label <- paste(ic99s$ic99, ' μg/mL', sep = '')

# creates plot for each line that overlays relative growth points and drm curve fit. also displays IC99 value.
pdf(file = paste('../figures/mic-', date,  '/rg+ic99_', cycle, '_bounded.pdf', sep = ''), height = 6, width = 12)
ggplot(data = rg_data, aes(x = as.numeric(abx_concentration), y = rg, color = isolate)) +
  geom_point(alpha = 0.5) +
  geom_line(linetype = 'dashed') +
  geom_line(data = predicted_rg, aes(x = as.numeric(abx_concentration), y = pred.growth, color = isolate)) +
  geom_text(data = ic99s, mapping = aes(x = 1, y = 5, label = ic99.label, group = isolate)) +
  geom_point(data = ic99s, shape = 4, size = 2, stroke = 2, mapping = aes(x = ic99, y = rg_at_IC99, group = isolate)) +
  scale_x_log10() +
  ggtitle(substitute(paste('Relative Growth and IC99 of WT and ', italic('fas Msm'), ' Strains Treated with Isoniazid'))) + xlab('isoniazid (μg/mL)') + ylab('relative growth') +
  ylim(0, 120) +
  plot_theme +
  scale_color_manual(values = colors) +
  theme(legend.position = 'none') +
  facet_wrap(~isolate, nrow = 2)
dev.off()


# creates data frame containing IC99 value for each isolate and writes to csv file.
IC_data <- data.frame(unlist(isolate_names), unlist(IC99_values))
names(IC_data)[1] <- "Isolate.Name"
names(IC_data)[2] <- "IC99.Value"

write_csv(IC_data, file = paste("../data/mic_", date, "/mic_values-", date, ".csv", sep = ""))
                   

