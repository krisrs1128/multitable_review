#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## An illustration of the geometric view of CCA using data from the WELL study.
##
## author: sankaran.kris@gmail.com
## date: 09/27/2017

###############################################################################
## Libraries and setup
###############################################################################
library("phyloseq")
library("tidyverse")
library("expm")
library("plot3D")

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
## Read and prepare data
###############################################################################
bc <- readRDS("../data/sample_data_bc.rds")
counts <- readRDS("../data/seqtab.rds")
colnames(counts) <- paste0("seq_", 1:ncol(counts))

bc_cols <- c("age", "Android_LM", "Total_FM")
bc_sub <- scale(bc[2:4, bc_cols])
counts_sub <- scale(asinh(counts[2:4, 1:3]))

x_grid <- seq(-2, 2, 0.1)
xyz_grid <- as.matrix(expand.grid(x_grid, x_grid, x_grid))
bc_plane <- matrix(nrow = nrow(xyz_grid), ncol = 3)
counts_plane <- matrix(nrow = nrow(xyz_grid), ncol = 3)
for (i in seq_len(nrow(xyz_grid))) {
  bc_plane[i, ] <- as.numeric(t(bc_sub) %*% xyz_grid[i, ])
  counts_plane[i, ] <- as.numeric(t(counts_sub) %*% xyz_grid[i, ])
}

sigma_bc <- cov(scale(bc[, bc_cols]))
sigma_counts <- cov(scale(counts[, 1:3]))

bc_circle <- matrix(nrow = 5000, ncol = 3)
counts_circle <- matrix(nrow = 5000, ncol = 3)
for (i in seq_len(5000)) {
  z <- rnorm(3)
  bc_circle[i, ] <- z / sqrt(t(z) %*% sigma_bc %*% z)
  counts_circle[i, ] <- z / sqrt(t(z) %*% sigma_counts %*% z)
}

## cca directions
sigma_cross <- cov(scale(bc[, bc_cols]))
cross_mat <- sqrtm(solve(sigma_bc)) %*% sigma_cross %*% sqrtm(solve(sigma_counts))
svd_cross <- svd(cross_mat)
cc1 <- cbind(
  "u" = sqrtm(solve(sigma_bc)) %*% svd_cross$u[, 1],
  "v" = sqrtm(solve(sigma_counts)) %*% svd_cross$v[, 1]
)

###############################################################################
## Plot the linear combinations of features
###############################################################################

png("../chapter/figure/cca/geometric.png")

## draw the cca points
arrows3D(
  theta = 40,
  phi = 20,
  ticktype = "detailed",
  xlab = rownames(bc_sub)[1],
  ylab = rownames(bc_sub)[2],
  zlab = rownames(bc_sub)[3],
  x0 = 0,
  y0 = 0,
  z0 = 0,
  x1 = cc1[1, 1],
  y1 = cc1[2, 1],
  z1 = cc1[3, 1],
  col = "#ff661a",
  xlim = c(-1, 1),
  ylim = c(-1, 1),
  zlim = c(-1, 1)
)

arrows3D(
  x0 = 0,
  y0 = 0,
  z0 = 0,
  x1 = cc1[1, 2],
  y1 = cc1[2, 2],
  z1 = cc1[3, 2],
  col = "#464282",
  add = TRUE
)

## draw the circle
points3D(
  x = bc_circle[, 1],
  y = bc_circle[, 2],
  z = bc_circle[, 3],
  col = "#ff9966",
  alpha = 0.6,
  cex = 0.2,
  add = TRUE
)

points3D(
  x = counts_circle[, 1],
  y = counts_circle[, 2],
  z = counts_circle[, 3],
  col = "#ccb3ff",
  alpha = 0.1,
  cex = 0.1,
  add = TRUE
)

## span of the body measurement data
points3D(
  x = bc_plane[, 1],
  y = bc_plane[, 2],
  z = bc_plane[, 3],
  col = "#ff9966",
  alpha = 0.1,
  cex = 0.1,
  add = TRUE
)

## span of the counts
points3D(
  x = counts_plane[, 1],
  y = counts_plane[, 2],
  z = counts_plane[, 3],
  col = "#ccb3ff",
  alpha = 0.1,
  cex = 0.1,
  add = TRUE
)

## labels for the body measurements
text3D(
  bc_sub[, 1],
  bc_sub[, 2],
  bc_sub[, 3],
  col = "#ff661a",
  labels = colnames(bc_sub),
  add = TRUE
)

## labels for the otu abundances
text3D(
  counts_sub[, 1],
  counts_sub[, 2],
  counts_sub[, 3],
  col = "#661aff",
  labels = colnames(counts_sub),
  add = TRUE
)
dev.off()
