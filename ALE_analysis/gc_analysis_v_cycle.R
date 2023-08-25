# 
# compiling growth curver data across cycles of MIC experiments
# 

library(tidyverse)
library(ggplot2)

colors <- c('#b486d2', '#a06dc1', '#844ca9', '#68308D', '#521b76', '#420d64', '#420d64',
            '#fbc068', '#f6ab3b', '#f4a732', '#F09814', '#f89500', '#d98200', '#d98200')


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

# collecting data about area under curve, carrying capacity, and growthr rate:
{
AUC_data1 <- read.csv("../data/mic_07262023/auc_gc_data-07262023.csv")
AUC_data2 <- read.csv("../data/mic_07282023/auc_gc_data-07282023.csv")
AUC_data3 <- read.csv("../data/mic_07312023/auc_gc_data-07312023.csv")
AUC_data4 <- read.csv("../data/mic_08022023/auc_gc_data-08022023.csv")
AUC_data5 <- read.csv("../data/mic_08042023/auc_gc_data-08042023.csv")
AUC_data6 <- read.csv("../data/mic_08062023/auc_gc_data-08062023.csv")
all_AUC_data <- do.call("rbind", list(AUC_data1, AUC_data2, AUC_data3, AUC_data4, AUC_data5, AUC_data6))
}

{
CC_data1 <- read.csv("../data/mic_07262023/carrycap_gc_data-07262023.csv")
CC_data2 <- read.csv("../data/mic_07282023/carrycap_gc_data-07282023.csv")
CC_data3 <- read.csv("../data/mic_07312023/carrycap_gc_data-07312023.csv")
CC_data4 <- read.csv("../data/mic_08022023/carrycap_gc_data-08022023.csv")
CC_data5 <- read.csv("../data/mic_08042023/carrycap_gc_data-08042023.csv")
CC_data6 <- read.csv("../data/mic_08062023/carrycap_gc_data-08062023.csv")
all_CC_data <- do.call("rbind", list(CC_data1, CC_data2, CC_data3, CC_data4, CC_data5, CC_data6))
}

{
R_data1 <- read.csv("../data/mic_07262023/rate_gc_data-07262023.csv")
R_data2 <- read.csv("../data/mic_07282023/rate_gc_data-07282023.csv")
R_data3 <- read.csv("../data/mic_07312023/rate_gc_data-07312023.csv")
R_data4 <- read.csv("../data/mic_08022023/rate_gc_data-08022023.csv")
R_data5 <- read.csv("../data/mic_08042023/rate_gc_data-08042023.csv")
R_data6 <- read.csv("../data/mic_08062023/rate_gc_data-08062023.csv")
all_R_data <- do.call("rbind", list(R_data1, R_data2, R_data3, R_data4, R_data5, R_data6))
}


# summarizing the data for each strain
gc_data_means_auc <- all_AUC_data %>%
  group_by(strain, abx_concentration, cycle) %>%
  summarize(average_auc = mean(average_auc), sd_auc = sd(average_auc))

gc_data_means_cc <- all_CC_data %>%
  group_by(strain, abx_concentration, cycle) %>%
  summarize(average_cc = mean(average_cc), sd_cc = sd(average_cc))

gc_data_means_rate <- all_R_data %>%
  group_by(strain, abx_concentration, cycle) %>%
  summarize(average_rate = mean(average_rate), sd_rate = sd(average_rate))

gc_data_means_auc$label <- as.factor(paste(gc_data_means_auc$strain, ', ', gc_data_means_auc$abx_concentration, ' µg/mL INH', sep = ''))
gc_data_means_auc$label <- factor(gc_data_means_auc$label, levels=c("fas, 0 µg/mL INH", "fas, 1 µg/mL INH", "fas, 2.5 µg/mL INH", "fas, 5 µg/mL INH", "fas, 10 µg/mL INH", "fas, 100 µg/mL INH", "fas, 1000 µg/mL INH", 
                                                                    "wt, 0 µg/mL INH", "wt, 1 µg/mL INH", "wt, 2.5 µg/mL INH", "wt, 5 µg/mL INH", "wt, 10 µg/mL INH", "wt, 100 µg/mL INH", "wt, 1000 µg/mL INH"))
gc_data_means_cc$label <- as.factor(paste(gc_data_means_cc$strain, ', ', gc_data_means_cc$abx_concentration, ' µg/mL INH', sep = ''))
gc_data_means_cc$label <- factor(gc_data_means_cc$label, levels=c("fas, 0 µg/mL INH", "fas, 1 µg/mL INH", "fas, 2.5 µg/mL INH", "fas, 5 µg/mL INH", "fas, 10 µg/mL INH", "fas, 100 µg/mL INH", "fas, 1000 µg/mL INH", 
                                                                    "wt, 0 µg/mL INH", "wt, 1 µg/mL INH", "wt, 2.5 µg/mL INH", "wt, 5 µg/mL INH", "wt, 10 µg/mL INH", "wt, 100 µg/mL INH", "wt, 1000 µg/mL INH"))
gc_data_means_rate$label <- as.factor(paste(gc_data_means_rate$strain, ', ', gc_data_means_rate$abx_concentration, ' µg/mL INH', sep = ''))
gc_data_means_rate$label <- factor(gc_data_means_rate$label, levels=c("fas, 0 µg/mL INH", "fas, 1 µg/mL INH", "fas, 2.5 µg/mL INH", "fas, 5 µg/mL INH", "fas, 10 µg/mL INH", "fas, 100 µg/mL INH", "fas, 1000 µg/mL INH", 
                                                                    "wt, 0 µg/mL INH", "wt, 1 µg/mL INH", "wt, 2.5 µg/mL INH", "wt, 5 µg/mL INH", "wt, 10 µg/mL INH", "wt, 100 µg/mL INH", "wt, 1000 µg/mL INH"))

# plots the compiled data into a line plot with the average auc, k, or r for each strain at each cycle
# cycle 1 had some inconsistent data due to slight differences in the protocol for that cycle. can use the commented out lines to remove cycle 1 from the plots.
pdf(file = '../figures/final_results/all_cycles_v_auc.pdf', height = 6, width = 12)
# ggplot(data = gc_data_means_auc %>% filter(cycle != 'C1'), aes(x = cycle, y = average_auc, group = label, color = label)) +
ggplot(data = gc_data_means_auc, aes(x = cycle, y = average_auc, group = label, color = label)) +
  ggtitle(substitute(paste('Average Area Under Growth Curve of WT and ', italic('fas Msm'), ' Strains Over Cycles of Isoniazid Treatment'))) + xlab('Cycle') + ylab('Average Area Under Curve') +
  geom_point(size=2) +
  geom_line(linewidth=1) +
  scale_color_manual(values = colors) +
  plot_theme +
  theme(legend.title = element_text(size = 0))
dev.off()

pdf(file = '../figures/final_results/all_cycles_v_carrycap.pdf', height = 6, width = 12)
# ggplot(data = gc_data_means_cc %>% filter(cycle != 'C1'), aes(x = cycle, y = average_cc, group = label, color = label)) +
ggplot(data = gc_data_means_cc, aes(x = cycle, y = average_cc, group = label, color = label)) +
    ggtitle(substitute(paste('Average Carrying Capacity of WT and ', italic('fas Msm'), ' Strains Over Cycles of Isoniazid Treatment'))) + xlab('Cycle') + ylab('Average Carrying Capacity') +
  geom_point(size=2) +
  geom_line(linewidth=1) +
  scale_color_manual(values = colors) +
  plot_theme +
  theme(legend.title = element_text(size = 0))
dev.off()

pdf(file = '../figures/final_results/all_cycles_v_rate.pdf', height = 6, width = 12)
# ggplot(data = gc_data_means_cc %>% filter(cycle != 'C1'), aes(x = cycle, y = average_cc, group = label, color = label)) +
ggplot(data = gc_data_means_rate, aes(x = cycle, y = average_rate, group = label, color = label)) +
  ggtitle(substitute(paste('Average Growth Rate of WT and ', italic('fas Msm'), ' Strains Over Cycles of Isoniazid Treatment'))) + xlab('Cycle') + ylab('Average Growth Rate') +
  geom_point(size=2) +
  geom_line(linewidth=1) +
  scale_color_manual(values = colors) +
  plot_theme +
  theme(legend.title = element_text(size = 0))
dev.off()
