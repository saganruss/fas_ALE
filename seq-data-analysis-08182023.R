# 
# script for analyzing reads per sample in sequencing data
# written by evan pepper
# 
{
  library(tidyverse)
  library(dplyr)
  library(ggplot2)
  }

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#################################################################################

plot.theme <- theme(title = element_text(size = 25),
                    axis.text = element_text(size = 16),
                    axis.title = element_text(size = 20),
                    legend.title = element_text(size = 20),
                    legend.text = element_text(size = 16))

# color vector for treatments
# treatment.colors <- c('purple', 'blue', 'red', 'gray')
#################################################################################
# reading in metadata files for my LLR experiments
# metadata <- read_csv(file = '../data/msm/08102023-metadata.csv')
# example.snps <- read_delim(file = '../data/msm/08182023-snps/DT1/snps.tab')
#################################################################################

sample.dirs <- list.dirs(path = '../data/snippy_results/08182023-ale-snps-0.05', recursive = F)
sample.ids <- unlist(lapply(strsplit(sample.dirs, '/'), '[', 5))

snp.table.list <- list()
for (d in 1:length(sample.dirs)) {
  sample.snps <- read_delim(file = paste(sample.dirs[d], '/snps.tab', sep = ''))
  sample.snps$sample.id <- sample.ids[d]
  snp.table.list[[d]] <- sample.snps
}

all.snps <- do.call(rbind, snp.table.list)

wt.snps <- all.snps %>% filter(sample.id == 'WT')

all.snps.filtered <- all.snps %>% filter(!POS %in% unique(wt.snps$POS))

write_csv(all.snps.filtered, file = paste("../data/snippy_results/snp_results/08182023_0.05.csv", sep = ""))
