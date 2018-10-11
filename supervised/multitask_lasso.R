#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## An experiment using lasso to model multiple responses in the WELL-China data.
##
## author: sankaran.kris@gmail.com
## date: 10/12/2017

###############################################################################
## Libraries and setup
###############################################################################
library("phyloseq")
library("tidyverse")
library("viridis")
library("glmnet")
source("../dimension_red/prep_tables.R")
source("../dimension_red/plot.R")

## cleaner ggplot theme
scale_colour_discrete <- function(...)
  scale_color_manual(
    values = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6', "#ffc100"),
    na.value = "#464646"
  )
scale_fill_discrete <- function(...)
  scale_fill_manual(
    values = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6',"#ffc100"),
    na.value = "#464646"
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
## Apply parallel lassos
###############################################################################
raw <- read_data()
opts <- list("filt_k" = 0.07)
processed <- process_data(
  raw$seqtab,
  raw$bc,
  raw$bc_full,
  raw$taxa,
  raw$tree
)
y <- sample_data(processed$ps) %>%
  select(-id, -number, -gender, -batch, -operator) %>%
  scale()
x <- scale(get_taxa(processed$ps))

n_lambda <- 30
lambdas <- seq(0.001, 0.7, length.out = n_lambda)
beta_hats <- array(dim = c(ncol(y), 1 + ncol(x), n_lambda))
fits <- list()

for (r in seq_len(ncol(y))) {
  message("tuning response ", r)
  fits[[r]] <- cv.glmnet(x, y[, r], lambda = lambdas, alpha = 0.3)
  beta_hats[r,, ] <- as.matrix(coef(fits[[r]]$glmnet.fit))
}

###############################################################################
## Plot the cross validation errors
###############################################################################
cv_err <- do.call(cbind, lapply(fits, function(x) x$cvm))
cv_err[cv_err > 1.4] <- NA
colnames(cv_err) <- colnames(y)
rownames(cv_err) <- lambdas

ggplot(melt(cv_err)) +
  geom_tile(
    aes(x = Var1, y = Var2, fill = value)
  ) +
  scale_fill_gradient2(midpoint = 1, low = "blue", high = "red") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0))

###############################################################################
## Plot the results
###############################################################################
seq_fam <- processed$mseqtab %>%
  select(seq_num, family) %>%
  unique()

species_order <- read_csv("species_order.csv") %>% unlist()
rownames(beta_hats) <- colnames(y)
colnames(beta_hats) <- c("intercept", colnames(x))
mbeta <- beta_hats %>%
  melt(
    varnames = c("feature", "seq_num", "lambda")
  ) %>%
  left_join(seq_fam) %>%
  mutate(
    feature = factor(feature, mass_ordering()),
    seq_num = factor(seq_num, c("intercept", species_order))
  ) %>%
  filter(seq_num != "intercept")

ggplot(mbeta) +
  geom_rect(
    aes(col = family),
    fill = "transparent", size = 0.2,
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) +
  geom_tile(
    aes(x = seq_num, y = lambda, fill = value),
    alpha = 0.7
  ) +
  scale_fill_gradient2(
    guide = guide_colorbar(ticks = FALSE, barheight = 0.5),
    mid = "#F8F8F8", low = "#40004b", high = "#00441b"
  ) +
  guides(
    col = guide_legend(
      override.aes = list(alpha = 1, size = 1.0),
      order = 1
    )
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_grid(feature ~ family, scale = "free", space = "free") +
  theme(
    panel.border = element_blank(),
    axis.text = element_blank(),
    panel.spacing = unit(0, "cm"),
    strip.text.y = element_text(hjust = 0, angle = 0),
    strip.text.x = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  "../chapter/figure/graph_lasso/multitask_lasso_hm_lambdas.png",
  width = 9.5,
  height = 6.5
)

## coefficient plot at just a single lambda
mbeta_sub <- mbeta %>%
  filter(lambda == 25)

ggplot(mbeta_sub) +
  geom_tile(
    aes(x = seq_num, y = feature, fill = value)
  ) +
  geom_rect(
    aes(col = family),
    fill = "transparent", size = 2,
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) +
  scale_colour_discrete() +
  scale_fill_gradient2(
    "Coef. ",
    guide = guide_colorbar(ticks = FALSE, barheight = 0.6),
    mid = "#F8F8F8", low = "#40004b", high = "#00441b"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  guides(color = guide_legend(nrow = 4, order = 1)) +
  facet_grid(. ~ family, scale = "free", space = "free") +
  theme(
    axis.text = element_blank(),
    panel.border = element_blank(),
    panel.spacing = unit(0, "cm"),
    axis.text.y = element_text(size = 6, angle = 0, hjust = 0),
    strip.text.x = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  "../chapter/figure/graph_lasso/multitask_lasso_hm.png",
  width = 5.9,
  height = 4.1
)
