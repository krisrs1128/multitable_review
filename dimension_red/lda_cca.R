#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## An experiment using Latent Dirichlet Allocation + Canonical Correlation
## Analysis simultaneously across count and gaussian tables.
##
## author: sankaran.kris@gmail.com
## date: 10/04/2017

###############################################################################
## Libraries and setup
###############################################################################
library("phyloseq")
library("rstan")
library("tidyverse")
library("reshape2")
library("ggrepel")
library("viridis")
library("PMA")
source("prep_tables.R")
source("plot.R")
source("lda_cca_plot.R")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

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
opts <- list(
  "vst" = FALSE,
  "stan_file" = "lda_cca.stan",
  "outdir" = "../chapter/figure/lda_cca/"
)
processed <- process_data(
  raw$seqtab,
  raw$bc,
  raw$bc_full,
  raw$taxa,
  raw$tree,
  opts
)

###############################################################################
## Run and plot LDA / CCA on the two (scaled) tables
###############################################################################
bc_mat <- processed$ps %>%
  sample_data() %>%
  select(-id, -number, -gender, -batch, -operator) %>%
  scale()

K <- 2
L1 <- 3
L2 <- 3

## initialize using results from sparse CCA
processed_vst <- process_data(raw$seqtab, raw$bc, raw$bc_full, raw$taxa, raw$tree)
cca_res <- CCA(
  scale(get_taxa(processed_vst$ps)),
  bc_mat,
  penaltyx = 0.6,
  penaltyz = 0.3,
  K = K
)

stan_data <- list(
  "n" = nrow(bc_mat),
  "p1" = ntaxa(processed$ps),
  "p2" = ncol(bc_mat),
  "K" = K,
  "L1" = L1,
  "L2" = L2,
  "sigma" = 1,
  "a0" = 1,
  "b0" = 1,
  "x" = get_taxa(processed$ps),
  "y" = bc_mat,
  "id_y" = diag(ncol(bc_mat)),
  "id_k" = diag(K),
  "id_l1" = diag(L1),
  "id_l2" = diag(L2),
  "zeros_k" = rep(0, K),
  "zeros_l1" = rep(0, L1),
  "zeros_l2" = rep(0, L2)
)

m <- stan_model(opts$stan_file)
init_vals <- list(
  "Bx" = cca_res$u,
  "By" = cca_res$v,
  "xi_s" = bc_mat %*% cca_res$v
)

vb_fit <- vb(m, data = stan_data, init = init_vals)
posterior <- rstan::extract(vb_fit)
rm(vb_fit)

###############################################################################
## Plot the results
###############################################################################
scv <- scale_color_viridis(
  guide = guide_colorbar(barheight = 0.15, ticks = FALSE)
)
seq_fam <- processed$mseqtab %>%
  seq_families()

mdist <- melt_parameters(posterior)
mdist$xi_s <- reshape_posterior_score(mdist$xi_s, raw$bc)
mdist$xi_x <- reshape_posterior_score(mdist$xi_x, raw$bc)
mdist$xi_y <- reshape_posterior_score(mdist$xi_y, raw$bc)
lda_cca_plots(mdist, seq_fam, processed, opts)
