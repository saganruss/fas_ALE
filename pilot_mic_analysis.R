# 
# analyzing MIC data from pilot analysis experiment started on 6/30/2023
# 

library(tidyverse)
library(ggplot2)

plot_theme <- theme(plot.title = element_text(hjust = 0.5),
                    panel.grid.major=element_line(colour="lightgray"),
                    panel.grid.minor=element_line(colour="lightgray"),
                    legend.key = element_rect(fill = "white"),
                    panel.background = element_rect(fill = "white")) 

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



annotated_data <- read_csv(file = "../data/20230703_ale_pilot/230630_annotateddata.CSV")

# modifying time to be in hours:
annotated_data$Time <- round(as.numeric(seconds(annotated_data$Time)) / 3600, digits = 3)

# plotting blank data:
blanks <-  startsWith(names(annotated_data), prefix = "blank")
blank_data <- cbind(annotated_data$Time, annotated_data[, blanks, drop = F])

blank_data <- blank_data %>% select(-c("blank...3"))

blank_data_plotting <- blank_data %>% pivot_longer(!names(blank_data[1]), names_to = "blank", values_to = "OD_600")
names(blank_data_plotting)[1] <- "Time"

ggplot(data = blank_data_plotting, aes(x = Time, y = OD_600, group = blank, color = blank)) +
  geom_point() +
  geom_line()

# plotting the sample data:
sample_data <- cbind(annotated_data[, !blanks, drop = F])

sample_data[, -1] <- sample_data[, -1] - rowMeans(blank_data[, -1])

{
sample_data_plotting <- sample_data %>% pivot_longer(!names(sample_data[1]), names_to = "sample", values_to = "OD_600")
sample_info <- strsplit(sample_data_plotting$sample, "_")
sample_data_plotting$strain <- unlist(lapply(sample_info, "[", 1))
sample_data_plotting$replicate <- unlist(lapply(sample_info, "[", 2))
sample_data_plotting$initial_OD <- as.numeric(unlist(lapply(sample_info, "[", 3)))
sample_data_plotting$abx_concentration <- as.numeric(unlist(lapply(sample_info, "[", 4)))
}

ggplot(data = sample_data_plotting, aes(x = Time, y = OD_600, group = sample, color = sample)) +
  geom_point() +
  geom_line()

sample_data_means <- sample_data_plotting %>% 
  group_by(Time, strain, initial_OD, abx_concentration) %>% 
  summarize(average_OD = mean(OD_600), max_OD = max(OD_600), min_OD = min(OD_600))

sample_data_means$sample_label <- paste(sample_data_means$strain, "_", sample_data_means$initial_OD, "_", sample_data_means$abx_concentration, sep = "")

curr_abx_conc <- 40
ggplot(data = sample_data_means %>% filter(abx_concentration == curr_abx_conc), aes(x = Time, y = average_OD, group = sample_label, color = sample_label)) +
  geom_point(size=2) +
  geom_line(size=1.3) +
  geom_errorbar(aes(ymax = max_OD, ymin = min_OD), size=0.8, alpha=0.6) +
  ggtitle(paste("Growth Curve of M. smegmatis Treated with ", curr_abx_conc, " μg/mL of INH", sep="")) +
  xlab("Time (hours)") + ylab(bquote(OD[600])) +
  plot_theme + 
  guides(color=guide_legend(bquote("Initial"~OD[600]))) +
  scale_colour_manual(labels = c("0.1", "0.25", "0.5", "", "", ""), 
                      values = c("royalblue", "blue", "darkblue", "royalblue", "blue", "darkblue")) +
  facet_wrap(. ~ strain)


od <- read_delim(file = '../data/od-measurements/od-07072023.txt', skip = 34, delim = '\t')
write_csv(od, file = '../data/od-measurements/od-07072023.csv')


