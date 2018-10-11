#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Illustration of canonical correlation analysis on the body composition and
## microbiome elements of the WELL-China data.
##
## author: sankaran.kris@gmail.com
## date: 09/26/2017

library("phyloseq")
library("tidyverse")
library("reshape2")
library("ggrepel")
library("viridis")
library("vegan")
source("prep_tables.R")
source("plot.R")

## cleaner ggplot theme
scale_colour_discrete <- function(...)
  scale_color_manual(
    values = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6', "#464646"),
    na.value = "black"
  )
scale_fill_discrete <- function(...)
  scale_fill_manual(
    values = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6', "#464646"),
    na.value = "black"
  )

theme_set(theme_bw())
theme_update(
  panel.background = element_rect(fill = "#F8F8F8"),
  panel.border = element_rect(size = 0.5),
  panel.grid = element_blank(),
  axis.ticks = element_blank(),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 6),
  axis.text = element_text(size = 6),
  axis.title = element_text(size = 8),
  strip.background = element_blank(),
  strip.text = element_text(size = 8),
  legend.key = element_blank()
)

###############################################################################
## read and prepare the data
###############################################################################
raw <- read_data()
opts <- list(filter_k = 0.5, filter_a = 5)
processed <- process_data(
  raw$seqtab,
  raw$bc,
  raw$bc_full,
  raw$taxa,
  raw$tree,
  opts
)

###############################################################################
## Run and plot CCA on the two (scaled) tables
###############################################################################
bc_mat <- processed$ps %>%
  sample_data() %>%
  select(-id, -number, -gender, -batch, -operator) %>%
  scale()
x_seq <- scale(get_taxa(processed$ps))
cca_res <- CCorA(bc_mat, x_seq)

## Plot the loadings
seq_fam <- processed$mseqtab %>%
  seq_families()

loadings <- prepare_loadings(
  list(cca_res$corr.Y.Cy, cca_res$corr.X.Cx),
  c("body_comp", "seq")
) %>%
  left_join(seq_fam)

plot_loadings(loadings, cca_res$Eigenvalues)
ggsave("../chapter/figure/cca/loadings.png", width = 4.4, height = 3.2)

## Plot the scores
scores <- prepare_scores(
  list(cca_res$Cx, cca_res$Cy),
  c("body_comp", "seq")
) %>%
  left_join(raw$bc)

mscores <- melt_scores(scores)
plot_scores(scores, "type", "Meas. Type", cca_res$Eigenvalues, c(-1, 1)) +
  link_scores(mscores) +
  scale_color_brewer(palette = "Set1")
ggsave("../chapter/figure/cca/scores_linked.png", width = 3.56, height = 2.6)

## color by weight
plot_scores(scores, "android_fm", "Android FM", cca_res$Eigenvalues, c(-1, 1)) +
  link_scores(mscores) +
  scale_color_viridis(
    guide = guide_colorbar(barheight = 3, barwidth = 0.15, ticks = FALSE)
  )
ggsave("../chapter/figure/cca/scores_android_fm.png", width = 4.56, height = 3.6)

## color by bacteroides / ruminoccocus ratios
scores <- scores %>%
  left_join(family_means(processed$mseqtab)) %>%
  mutate(rl_ratio = tanh(rl_ratio))
plot_scores(scores, "rl_ratio", "tanh(Bact. - Rumino.)", cca_res$Eigenvalues) +
  link_scores(mscores) +
  scale_color_viridis(
    guide = guide_colorbar(barheight = 3, barwidth = 0.15, ticks = FALSE)
  )
ggsave("../chapter/figure/cca/scores_rl_ratio.png", width = 3.67, height = 3.4)
