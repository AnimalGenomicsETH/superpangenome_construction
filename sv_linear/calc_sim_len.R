#!/usr/bin/env Rscript
# loading library and ggplot preference
library(tidyverse)
library(patchwork)

theme_set(theme_bw(base_size = 24) +
        theme(panel.grid = element_blank()))
colpal <- c(
        "#E69F00", "#56B4E9", "#009E73", "#D55E00",
        "#CC79A7", "#0072B2", "#F0E442", "#999999"
)

options(
        ggplot2.discrete.colour = colpal,
        ggplot2.discrete.fill = colpal,
        ggplot2.continuous.colour = colpal,
        ggplot2.continuous.fill = colpal
)

# set working directory
setwd("/cluster/work/pausch/danang/psd/scratch/sim/lr_sim/simulated_lr")

datont <- read.table("ont_len.tsv", header = FALSE, stringsAsFactors = FALSE)
datpb <- read.table("pacbio_len.tsv", header = FALSE, stringsAsFactors = FALSE)
datall <- rbind(datont, datpb)
head(datall)
colnames(datall) <- c("prog", "rep", "noread", "readlen")
table(datall$rep)

datall %>%
        group_by(prog, rep) %>%
        summarise(
                totread = n(),
                mlen = mean(readlen),
                maxlen = max(readlen),
                minlen = min(readlen)
        )


pl1 <- ggplot(datall, aes(x = readlen, fill = prog)) +
        geom_density(alpha = 0.4) +
        theme(legend.position = "bottom")
pl1

ggsave("read_dist.png")