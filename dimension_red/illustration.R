#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Just some pedagogical examples for the survey.
##
## author: sankaran.kris@gmail.com
## date: 09/24/2017

###############################################################################
## libraries
###############################################################################
library("tidyverse")
scale_colour_discrete <- function(...)
  scale_colour_brewer(..., palette="Set2")
scale_fill_discrete <- function(...)
  scale_fill_brewer(..., palette="Set2")

theme_set(theme_bw())
theme_update(
  panel.border = element_rect(size = 0.7),
  panel.grid = element_blank(),
  axis.ticks = element_blank(),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 6),
  axis.text = element_text(size = 8),
  axis.title = element_text(size = 10),
  strip.background = element_blank(),
  strip.text = element_text(size = 8),
  legend.key = element_blank()
)

project <- function(x, pc_x, K = 2) {
  v <- pc_x$rotation
  x %*% v[, 1:K] %*% t(v[, 1:K])
}

proj_plot <- function(bc, proj_x) {
  ggplot() +
    geom_segment(
      data = data.frame("proj" = proj_x, bc),
      aes(x = weight_dxa, xend = proj.weight_dxa, y = Total_FM, yend = proj.Total_FM),
      alpha = 0.6, col = "#9370db"
    ) +
    geom_point(
      data = data.frame(bc),
      aes(x = weight_dxa, y = Total_FM),
      alpha = 0.8
    ) +
    geom_point(
      data = data.frame(proj_x),
      aes(x = weight_dxa, y = Total_FM),
      col = "#9370db", alpha = 0.3
    ) +
    coord_fixed() +
    labs(x = "Weight", y = "Fat Mass")
}

###############################################################################
## run pc examples
###############################################################################
bc <- readRDS("../data/sample_data_bc.rds") %>%
  filter(gender == "Female") %>%
  select(-id, -Number, -gender) %>%
  as.matrix() %>%
  scale()

bc_sub <- bc[, c("weight_dxa", "Total_FM")]
pc_x <- prcomp(bc_sub)
proj_x <- project(bc_sub, pc_x, 1)
proj_plot(bc_sub, proj_x)
ggsave("../chapter/figure/pca/proj_plot_1.png", width = 3, height = 3)

pc_x <- prcomp(bc)
proj_x <- project(bc, pc_x, 2)
proj_plot(bc, proj_x)
ggsave("../chapter/figure/pca/proj_plot_2.png", width = 3, height = 3)

pc_x <- prcomp(bc)
proj_x <- project(bc, pc_x, 1)
proj_plot(bc, proj_x)
ggsave("../chapter/figure/pca/proj_plot_3.png", width = 3, height = 3)

###############################################################################
## maximum variance interpretation
###############################################################################
z <- rnorm(ncol(bc))
rand_proj <- bc %*% z
rand_proj <- rand_proj / sqrt(sum(rand_proj ^ 2))
rand_proj <- as.numeric(rand_proj)

contrasts <- data_frame(
  i = 1:nrow(bc),
  means = rowMeans(bc),
  pc1 = pc_x$x[, "PC1"],
  pc2 = pc_x$x[, "PC2"],
  random = rand_proj
) %>%
  gather(var, value, -i) %>%
  mutate(
    var = factor(var, levels = c("pc1", "pc2", "means", "random"))
  )

ggplot(contrasts) +
  geom_histogram(
    aes(x = value),
    bins = 70
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_grid(var ~ .)
ggsave("../chapter/figure/pca/var_plot.png", width = 4, height = 3)
