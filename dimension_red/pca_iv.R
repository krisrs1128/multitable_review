#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Illustration of PCA-IV on the body composition and microbiome elements of the
## WELL-China data.
##
## author: sankaran.kris@gmail.com
## date: 09/29/2017

library("phyloseq")
library("tidyverse")
library("reshape2")
library("ggrepel")
library("viridis")
library("ade4")
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
## Run PCA-IV
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

bc_mat <- processed$ps %>%
  sample_data() %>%
  select(-id, -number, -gender, -batch, -operator) %>%
  scale()

K <- 3
pca_micro <- dudi.pca(
  scale(get_taxa(processed$ps)),
  center = FALSE,
  scale = FALSE,
  scan = FALSE,
  nf = K
)
pcaiv_res <- pcaiv(
  pca_micro, bc_mat,
  scan = FALSE,
  nf = K
)

###############################################################################
## study equivalent of scores (projections of x and y onto principal axes)
###############################################################################
pcaiv_res
summary(pcaiv_res)
plot(pcaiv_res)

## correlate columns from both data frames with the principal axes
loadings <- prepare_loadings(
  list(0.01 * pcaiv_res$fa, pcaiv_res$c1),
    c("body_comp", "seq")
) %>%
  left_join(processed$mseqtab)
plot_loadings(loadings, pcaiv_res$eig, size_breaks=c(-6, 6)) +
  scale_size_continuous(range = c(0.4, 3), breaks=c(-3, 3))
ggsave("../chapter/figure/pca_iv/loadings.png", width = 4.6, height = 3.0)

## project the samples onto the principal axes
scores <- prepare_scores(
  list(pcaiv_res$li, pcaiv_res$ls),
  c("body_comp", "seq")
) %>%
  left_join(raw$bc)

mscores <- melt_scores(scores)
plot_scores(scores, "android_fm", "Android FM", pcaiv_res$eig) +
  link_scores(mscores) +
  scale_color_viridis(
    guide = guide_colorbar(barheight = 0.15, ticks = FALSE)
  )
ggsave("../chapter/figure/pca_iv/scores_android_fm.png", width = 3.56, height = 1.8)

scores <- scores %>%
  left_join(family_means(processed$mseqtab))
plot_scores(scores, "rl_ratio", "tanh(Bact. - Rumino.)", pcaiv_res$eig) +
  link_scores(mscores) +
  scale_color_viridis(
    guide = guide_colorbar(barheight = 0.15, ticks = FALSE)
  )
ggsave("../chapter/figure/pca_iv/scores_rl_ratio.png", width = 3.56, height = 1.8)
