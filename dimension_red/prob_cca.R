#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Applying probabilistic CCA to real data.
##
## author: sankaran.kris@gmail.com
## date: 10/04/2017

###############################################################################
## Libraries and setup
###############################################################################
library("tidyverse")
library("rstan")
library("viridis")
source("prep_tables.R")
source("plot.R")
set.seed(9282017)

scale_colour_discrete <- function(...)
  scale_colour_brewer(..., palette="Set2")
scale_fill_discrete <- function(...)
  scale_fill_brewer(..., palette="Set2")

theme_set(theme_bw())
theme_update(
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
## Real data application
###############################################################################
## read and prepare data. Note that it would be possible to use more taxa
## (unlike in ordinary CCA), but this gets much slower. A faster alternative
## that does use all the species is in lda_cca.R
raw <- read_data()
processed <- process_data(raw$seqtab, raw$bc, raw$taxa)

K <- 2
L1 <- 3
L2 <- 3
stan_data <- list(
  "K" = K,
  "L1" = L1,
  "L2" = L2,
  "p1" = ncol(processed$bc),
  "p2" = ncol(processed$x_seq),
  "n" = nrow(processed$bc),
  "tau" = 5,
  "sigma_x" = 1,
  "sigma_y" = 1,
  "x" = scale(processed$bc),
  "y" = scale(processed$x_seq),
  "id_x" = diag(ncol(processed$bc)),
  "id_y" = diag(ncol(processed$x_seq)),
  "id_k" = diag(K),
  "id_l1" = diag(L1),
  "id_l2" = diag(L2),
  "zeros_k" = rep(0, K),
  "zeros_l1" = rep(0, L1),
  "zeros_l2" = rep(0, L2)
)

m <- stan_model("prob_cca.stan")
vb_fit <- vb(m, stan_data)
posterior <- rstan::extract(vb_fit)
rm(vb_fit)

pmean <- parameter_means(posterior)
shared_scores <- cbind(pmean$xi_s, processed$bc)
colnames(shared_scores)[1:K] <- paste0("Axis", 1:K)

## should try to plot posteriors, not just means
ggplot(shared_scores) +
  geom_point(
    aes(x = Axis1, y = Axis2, col = weight_dxa)
  ) +
  scale_color_viridis(
    "Weight ",
    guide = guide_colorbar(barheight= 0.15, ticks = FALSE)
  )

rownames(pmean$Wx) <- colnames(processed$bc)
rownames(pmean$Wy) <- colnames(processed$x_seq)
rownames(pmean$Bx) <- colnames(processed$bc)
rownames(pmean$By) <- colnames(processed$x_seq)

seq_families <- processed$mseqtab %>%
  select(seq_num, family) %>%
  unique()

loadings_within <- prepare_loadings(
  list(pmean$Wx, pmean$Wy),
  c("body_comp", "seq")
) %>%
  left_join(seq_families)
plot_loadings(loadings_within, c(1, 1))

loadings_between <- prepare_loadings(
  list(cbind(pmean$Bx, 1), cbind(pmean$By, 1)),
  c("body_comp", "seq"),
) %>%
  left_join(seq_families)
plot_loadings(loadings_between, c(1, 1)) +
  scale_size(range = c(0, 4), guide = FALSE)
