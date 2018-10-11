#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## An experiment using Latent Dirichlet Allocation + Canonical Correlation
## Analysis simultaneously across count and gaussian tables. The only difference
## with lda_cca.R is that this uses a dirichlet prior for the scores.
##
## author: sankaran.kris@gmail.com
## date: 10/05/2017

###############################################################################
## Libraries and setup
###############################################################################
library("phyloseq")
library("rstan")
library("tidyverse")
library("reshape2")
library("ggrepel")
library("viridis")
source("prep_tables.R")
source("plot.R")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## cleaner ggplot theme
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
## read and prepare the data
###############################################################################
raw <- read_data()
opts <- list("filt_k" = 0.02, "rlog" = FALSE)
processed <- process_data(raw$seqtab, raw$bc, raw$taxa, opts)

###############################################################################
## Run and plot LDA / CCA on the two (scaled) tables
###############################################################################
bc_mat <- scale(processed$bc)

K <- 2
L1 <- 4
L2 <- 2

stan_data <- list(
  "n" = nrow(bc_mat),
  "p1" = ncol(processed$x_seq),
  "p2" = ncol(bc_mat),
  "K" = K,
  "L1" = L1,
  "L2" = L2,
  "sigma" = 1,
  "id_y" = diag(ncol(bc_mat)),
  "alpha_k" = rep(1, K),
  "alpha_l1" = rep(1, L1),
  "alpha_l2" = rep(1, L2),
  "gamma" = rep(1, ncol(processed$x_seq)),
  "x" = processed$x_seq,
  "y" = bc_mat
)

m <- stan_model("lda_cca_dir.stan")
vb_fit <- vb(m, data = stan_data)

posterior <- rstan::extract(vb_fit)
pmeans <- parameter_means(posterior)
for (i in seq_along(pmeans)) {
  pmeans[[i]] <- cbind(pmeans[[i]], 1) # in case only used 2 dimensions
}

###############################################################################
## Plot the results
###############################################################################
scv <- scale_color_viridis(
  guide = guide_colorbar(barheight = 0.15, ticks = FALSE)
)
seq_families <- processed$mseqtab %>%
  select(seq_num, family) %>%
  unique()

## Plot some of the shared scores and loadings
plot_scores_wrapper(pmeans$xi_s, raw, processed, scv)

rownames(pmeans$By) <- colnames(processed$bc)
rownames(pmeans$Bx) <- colnames(processed$x_seq)
loadings <- prepare_loadings(
  list(pmeans$By, pmeans$Bx),
  c("body_comp", "seq")
) %>%
  left_join(seq_families)

plot_loadings(
  loadings %>% filter(type == "seq", !is.na(family)),
  c(1, 1),
  a = 0.4
) +
  facet_wrap(~family, ncol = 4) +
  scale_color_brewer(guide = FALSE) +
  scale_size_continuous(range = c(0, 2), guide = FALSE)
plot_loadings(loadings %>% filter(type == "body_comp"), c(1, 1))

## now plot unshared scores (species abundances first)
plot_scores_wrapper(pmeans$xi_x, raw, processed, scv)
plot_scores_wrapper(pmeans$xi_y, raw, processed, scv)

## Now plot unshared loadings
rownames(pmeans$Wx) <- colnames(processed$x_seq)
loadings_x <- prepare_loadings(list(pmeans$Wx), "seq") %>%
  left_join(seq_families) %>%
  filter(!is.na(family))
plot_loadings(loadings_x, c(1, 1), a = 0.4) +
  facet_wrap(~family, ncol = 4) +
  scale_color_brewer(palette = "Set2", guide = FALSE) +
  scale_size_continuous(range = c(0, 2), guide = FALSE)

rownames(pmeans$Wy) <- colnames(processed$bc)
loadings_y <- prepare_loadings(list(pmeans$Wy), "body_comp") %>%
  mutate(seq_num = "NA") %>%
  left_join(seq_families)
plot_loadings(loadings_y, c(1, 1)) +
  scale_size_continuous(range = c(1, 4), guide = FALSE)
