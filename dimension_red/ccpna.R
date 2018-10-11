#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Application of the canonical correspondence analysis to the body composition
## data.
##
## author: sankaran.kris@gmail.com
## date: 10/06/2017

###############################################################################
## Libraries and setup
###############################################################################
library("phyloseq")
library("tidyverse")
library("vegan")
library("viridis")
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
## Apply canonical correspondence analysis
###############################################################################
raw <- read_data()
opts <- list(vst = FALSE)
processed <- process_data(
  raw$seqtab,
  raw$bc,
  raw$bc_full,
  raw$taxa,
  raw$tree,
  opts
)
x <- processed$ps %>%
  sample_data() %>%
  select(-id, -number, -gender, -batch, -operator)
y <- get_taxa(processed$ps)
cca_res <- cca(
  y ~ age + weight_dxa + android_fm + gynoid_fm + aoi,
  data.frame(x)
)

###############################################################################
## Plot loadings and scores
###############################################################################
seq_fam <- processed$mseqtab %>%
  seq_families()

## Plot the scores
cc_scores <- prepare_scores(
  scores(cca_res, choices = 1:3)[2],
  "body_comp"
) %>%
  left_join(raw$bc)

plot_scores(cc_scores, "android_fm", "Android FM", cca_res$CCA$eig) +
  scale_color_viridis(
    guide = guide_colorbar(barheight = 0.15, ticks = FALSE)
  )
ggsave("../chapter/figure/ccpna/scores_total_lm.png", width = 5.47, height = 3.19)

loadings <- prepare_loadings(
  list(4 * cca_res$CCA$biplot, cca_res$CCA$v),
  c("body_comp", "seq")
) %>%
  left_join(seq_fam)
plot_loadings(loadings, cca_res$CCA$eig, c(-8, 4)) +
  scale_size_continuous(range = c(0.1, 3))
ggsave("../chapter/figure/ccpna/loadings.png", width = 9.81, height = 4.62)
