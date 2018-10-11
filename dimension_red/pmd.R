#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Application of the penalized matrix decomposition to the body composition
## data.
##
## author: sankaran.kris@gmail.com
## date: 10/06/2017

###############################################################################
## Libraries and setup
###############################################################################
library("phyloseq")
library("tidyverse")
library("PMA")
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
## read and prepare the data
###############################################################################
raw <- read_data()
processed <- process_data(raw$seqtab, raw$bc, raw$bc_full, raw$taxa, raw$tree)
x <- sample_data(processed$ps) %>%
  select_if(is.numeric) %>%
  scale()
y <- scale(get_taxa(processed$ps))
cca_res <- CCA(x, y, penaltyx = 0.7, penaltyz = 0.3, K = 3)

## Plot the scores
scores <- prepare_scores(
  list(x %*% cca_res$u, y %*% cca_res$v),
  c("body_comp", "seq")
) %>%
  left_join(sample_data(processed$ps))

mscores <- melt_scores(scores)
plot_scores(scores, "type", "Meas. Type", cca_res$d) +
  link_scores(mscores) +
  scale_color_brewer(palette = "Set1")

## color by android fat mass
plot_scores(scores, "android_fm", "android fm", cca_res$d, c(-3, 3)) +
  link_scores(mscores) +
  scale_color_viridis(
    guide = guide_colorbar(barheight = 3, barwidth = 0.15, ticks = FALSE)
  )
ggsave("../chapter/figure/pmd/scores_android_fm.png", width = 4.4, height = 3.1)

## color by ruminoccocus / lachnospiraceae ratios
scores <- scores %>%
  left_join(family_means(processed$mseqtab))
plot_scores(scores, "rl_ratio", "tanh(Bact. - Rumino.)", cca_res$d) +
  link_scores(mscores) +
  scale_color_viridis(
    guide = guide_colorbar(barheight = 3, barwidth = 0.15, ticks = FALSE)
  )

## Plot the loadings
rownames(cca_res$u) <- colnames(x)
rownames(cca_res$v) <- colnames(y)
seq_fam <- processed$mseqtab %>%
  seq_families()

loadings <- prepare_loadings(
  list(data.frame(cca_res$u), data.frame(cca_res$v)),
  c("body_comp", "seq")
) %>%
  left_join(seq_fam)

large_species <- loadings %>%
  filter(
    sqrt(Axis.1 ^ 2 + Axis.2 ^ 2) > 0.15,
    type == "seq"
  )

plot_loadings(
  loadings,
  cca_res$d
) +
  geom_text_repel(
    data = large_species,
    aes(
      x = Axis.1,
      y = Axis.2,
      label = seq_num,
      col = family
    ),
    size = 1.5
  ) +
ggsave("../chapter/figure/pmd/loadings.png", width = 4.5, height = 3.3)

mlarge_species <- melt(
  data.frame(
    "number" = rownames(x),
    y[, large_species$seq_num],
    x[, c("android_fm", "gynoid_fm")]
  ),
  id.vars = c("number", "android_fm", "gynoid_fm"),
  variable.name = "seq_num"
) %>%
  left_join(seq_fam)

mlarge_species$seq_num <- factor(
  mlarge_species$seq_num,
  loadings %>%
  arrange(dplyr::desc(Axis.1)) %>%
    .[["seq_num"]]
)

ggplot(mlarge_species) +
  geom_hline(yintercept = 0, size = 0.1, alpha = 0.8) +
  geom_vline(xintercept = 0, size = 0.1, alpha = 0.8) +
  stat_smooth(
    aes(x = value, y = android_fm, col = family),
    fill = "#adadad",
    size = 0.5,
    method = "lm"
  ) +
  geom_point(
    aes(x = value, y = android_fm, col = family),
    size = 0.7, alpha = 0.8
  ) +
  facet_wrap(~seq_num, ncol = 8) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  "../chapter/figure/pmd/android_fm_species.png",
  width = 7.06,
  height = 4.41
)
