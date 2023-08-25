# 
# compiling MIC data across cycles of MIC experiments
# 

library(tidyverse)
library(ggplot2)

colors <- c('#b486d2', '#a06dc1', '#844ca9', '#68308D', '#521b76', '#420d64',
            '#fbc068', '#f6ab3b', '#f4a732', '#F09814', '#f89500', '#d98200')


plot_theme <- theme(plot.title = element_text(hjust = 0.5),
                    panel.grid.major=element_line(colour="lightgray"),
                    panel.grid.minor=element_line(colour="lightgray"),
                    legend.key = element_rect(fill = "white"),
                    panel.background = element_rect(fill = "white"),
                    axis.text.x = element_text(size = 15),
                    axis.text.y = element_text(size = 15),
                    legend.text = element_text(size = 15),
                    legend.title = element_text(size = 15),
                    axis.title = element_text(size = 15))

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# collecting data from first iteration of experiment:
# IC_data1 <- read.csv("../data/mic_07072023/mic_values-07072023_unbounded.csv")
# IC_data2 <- read.csv("../data/mic_07112023/mic_values-07112023_unbounded.csv")
# IC_data3 <- read.csv("../data/mic_07142023/mic_values-07142023_unbounded.csv")
# IC_data4 <- read.csv("../data/mic_07182023/mic_values-07182023_unbounded.csv")
# IC_data5 <- read.csv("../data/mic_07212023/mic_values-07212023_unbounded.csv")
# all_data <- do.call("rbind", list(IC_data1, IC_data2, IC_data3, IC_data4, IC_data5))


# collecting data from second iteration of experiment:
IC_data1 <- read.csv("../data/mic_07262023/mic_values-07262023.csv")
IC_data2 <- read.csv("../data/mic_07282023/mic_values-07282023.csv")
IC_data3 <- read.csv("../data/mic_07312023/mic_values-07312023.csv")
IC_data4 <- read.csv("../data/mic_08022023/mic_values-08022023.csv")
IC_data5 <- read.csv("../data/mic_08042023/mic_values-08042023.csv")
IC_data6 <- read.csv("../data/mic_08062023/mic_values-08062023.csv")
all_data <- do.call("rbind", list(IC_data1, IC_data2, IC_data3, IC_data4, IC_data5, IC_data6))

{
data_info <- strsplit(all_data$Isolate.Name, "_")
all_data$strain <- as.factor(unlist(lapply(data_info, "[", 1)))
all_data$line <- unlist(lapply(data_info, "[", 2))
all_data$cycle <- unlist(lapply(lapply(data_info, "[", 3), substring, 2))
all_data$Line <- paste(all_data$strain, "_", all_data$line, sep = "")
all_data$foldchange <- 0
}

# calculates fold change (Note: change 'C0' to 'C1' when processing data from first iteration)
for (i in 0:nrow(all_data)) {
  line <- all_data[i, ] %>% pull(Line)
  initial_line <- paste(line, "_C1", sep = "")
  initial_IC99 <- all_data %>% filter(Isolate.Name == initial_line) %>% pull(IC99.Value)
  curr_IC99 <- all_data[i,] %>% pull(IC99.Value)
  
  all_data[i, ]$foldchange <- curr_IC99 / initial_IC99
  
}

# adding another column for plotting label
all_data$label <- paste(all_data$strain, ', line ', substr(all_data$line, 2, 3), sep = '')

# plots the ic99 value at each cycle for the entire ALE experiment
pdf(file = '../figures/final_results/all_cyclesvic99.pdf', height = 6, width = 12)
ggplot(data = all_data, aes(x = cycle, y = IC99.Value, group = label, color = label)) +
  ggtitle(substitute(paste('IC99 of WT and ', italic('fas Msm'), ' Strains Over Cycles of Isoniazid Treatment'))) + xlab('Cycle') + ylab('IC99 (ug/mL)') +
  geom_point(size=2) +
  geom_line(linewidth=1) +
  scale_color_manual(values = colors) +
  plot_theme +
  theme(legend.title = element_text(size = 0))
dev.off()

# plots the FC of the ic99 value at each cycle compared to the initial ic99 value for the entire ALE experiment
pdf(file = '../figures/mic-07212023//all_cyclesvic99FC.pdf', height = 6, width = 12)
ggplot(data = all_data, aes(x = cycle, y = foldchange, group = label, color = label)) +
  ggtitle(substitute(paste('MIC Fold Change of WT and ', italic('fas Msm'), ' Strains Over Cycles of Isoniazid Treatment'))) + xlab('Cycle') + ylab('MIC Fold Change') +
  geom_point(size=2) +
  geom_line(linewidth=1) +
  scale_color_manual(values = colors) +
  plot_theme +
  theme(legend.title = element_text(size = 0))
dev.off()


# creates array with just the samples that were sequenced
sequenced_samples <- c(1, 25, 7, 31, 2, 3, 5, 6, 8, 9, 11, 12, 14, 15, 17, 18, 20, 21, 23, 24, 26, 27, 29, 30, 32, 33, 35, 36)
sequenced_data <- all_data
sequenced_data <- all_data[sequenced_samples, ]

# creates array with just the fas samples that were sequenced
fas_sequenced_samples <- c(7, 31, 8, 9, 11, 12, 20, 21, 23, 24, 32, 33, 35, 36)
fas_sequenced_data <- all_data[fas_sequenced_samples, ]

# creates array with just the WT samples that were sequenced
wt_sequenced_samples <- c(1, 25, 2, 3, 5, 6, 14, 15, 17, 18, 26, 27, 29, 30)
wt_sequenced_data <- all_data[wt_sequenced_samples, ]

colors <- c('#b486d2', '#a06dc1', '#844ca9', '#68308D', '#521b76',
            '#fbc068', '#f6ab3b', '#f4a732', '#F09814', '#f89500')

# plots the fold change of the ic99 of fas sequenced samples at each cycle compared to their initial ic99 value
pdf(file = '../figures/final_results/fas_sequenced_all_cyclesvic99FC.pdf', height = 6, width = 12)
ggplot(data = fas_sequenced_data, aes(x = cycle, y = foldchange, group = label, color = label)) +
  ggtitle(substitute(paste('MIC Fold Change of Sequenced ', italic('fas Msm'), ' Strains Over Cycles of Isoniazid Treatment'))) + xlab('Cycle') + ylab('MIC Fold Change') +
  geom_point(size=2) +
  geom_line(linewidth=1) +
  scale_color_manual(values = colors) +
  plot_theme +
  theme(legend.title = element_text(size = 0))
dev.off()


