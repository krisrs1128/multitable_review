
################################################################################
# Simulate data according to Daniela Witten's PMD paper
################################################################################

###############################################################################
## Libraries and setup
###############################################################################

library("reshape2")
library("tidyverse")
library("PMA")
set.seed(04032016)

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

#' Model with common sources, but different scores for each
#'
#' @param W A list of matrices containing scores for each table.
#' @param S A single source matrix underlying all the tables.
#' @return A list of matrices, each with the same sources S but different
#' scores W.
common_source_model <- function(W, S, opts) {
  X <- list()
  n <- nrow(W[[1]])
  p <- nrow(S)

  for(l in seq_along(W)) {
    X[[l]] <- W[[l]] %*% t(S) + matnorm(n, p, opts$sigma)
  }

  X
}

#' i.i.d. Gaussian Matrix
matnorm <- function(n, p, mu = 0, sigma = 1) {
  matrix(
    rnorm(n * p, mu, sigma),
    n, p
  )
}

norm <- function(x) {
  sqrt(sum(x ^ 2))
}

#' Plot simulation experiment
plot_simulation <- function(mW, pmd_res) {
  mW_hat <- melt(pmd_res$ws)
  colnames(mW_hat) <- c("i", "k", "w", "table")
  mW_hat$table[mW_hat$table == 1] <- "X"
  mW_hat$table[mW_hat$table == 2] <- "Y"
  mW_hat$type <- "recovered"
  mW_hat$k <- as.factor(mW_hat$k)
  mW$k <- as.factor(mW$k)

  plots <- list()
  plots[["sequence"]] <- ggplot() +
    geom_line(
      data = mW,
      aes(x = i, y = w, col = k)
    ) +
    geom_point(
      data = mW_hat,
      aes(x = i, y = w, col = k),
      alpha = 0.6,
      size = 1
    ) +
    facet_grid(table ~ .)

  mW_bind <- bind_rows(mW, mW_hat)
  mW_bind$w <- mW_bind$w + runif(nrow(mW_bind), -0.01, 0.01)
  mW_points <- mW_bind %>%
    spread(k, w) %>%
    dplyr::rename(V1 = `1`, V2 = `2`)
  mW_links <- mW_bind %>%
    dcast(i + table ~ k + type, value.var = "w")

  plots[["scatter"]] <- ggplot(mW_bind) +
    geom_hline(yintercept = 0, alpha = 0.2) +
    geom_vline(xintercept = 0, alpha = 0.2) +
    geom_point(
      data = mW_points,
      aes(x = V1, y = V2, col = type, shape = table),
      size = 1,
      alpha = 0.6
    ) +
    geom_segment(
      data = mW_links,
      aes(x = `1_truth`, xend = `1_recovered`, y = `2_truth`, yend = `2_recovered`),
      size = 0.2,
      alpha = 0.1
    ) +
    scale_color_brewer(
      palette = "Set2",
      guide = guide_legend(override.aes = list(alpha = 1, size = 2))
    )

  plots
}

opts <- list(
  "n" = 504, # want divisible by 8
  "p" = 20,
  "k" = 2,
  "l" = 2,
  "sigma0" = 2,
  "sigma" = 0.01
)

###############################################################################
## Simulate data
###############################################################################
S <- 500 * qr.Q(qr(matnorm(opts$p, opts$k, opts$sigma0))) # same source between the two matrices
W <- replicate(opts$l, matrix(0, opts$n, opts$k), simplify = F)
W[[1]][, 1] <- c(
  rep(0, opts$n / 8),
  rep(-10, opts$n / 8),
  rep(10, opts$n / 8),
  rep(-10, opts$n / 8),
  rep(0, opts$n / 2)
)
W[[1]][, 2] <- c(
  rep(10, opts$n / 4),
  rep(-10, opts$n / 4),
  rep(0, opts$n / 2)
)
W[[2]][, 1] <- c(
  rep(0, opts$n / 4),
  rep(10, opts$n / 2),
  rep(-10, opts$n / 4)
)
W[[2]][, 2] <- c(
  rep(-10, opts$n / 4),
  rep(0, opts$n / 4),
  rep(10, opts$n / 4),
  rep(0, opts$n / 4)
)

for (i in seq_along(W)) {
  for (j in seq_len(ncol(W[[i]]))) {
    W[[i]][, j] <- W[[i]][, j] / norm(W[[i]][, j])
  }
}

###############################################################################
## Plot the simulated data
###############################################################################
mW <- melt(W)
colnames(mW) <- c("i", "k", "w", "table")
mW$table[mW$table == 1] <- "X"
mW$table[mW$table == 2] <- "Y"
mW$type <- "truth"

###############################################################################
## Simulation according to common source model
###############################################################################
X <- common_source_model(W, S, opts)
colnames(X[[2]]) <- paste0("Y", seq_len(ncol(X[[2]])))

x_ix <- sample(seq_len(ncol(X[[1]])), 4)
y_ix <- sample(seq_len(ncol(X[[1]])), 4)
pairs(X[[1]][, x_ix], asp = 1, main = "Four columns of X")
pairs(X[[2]][, y_ix], asp = 1, main = "Four columns of Y")

pairs(cbind(X[[1]][, x_ix[1:2]], X[[2]][, y_ix[1:2]]), asp = 1,
      main = "Two columns of X vs. Two columns of Y")

pmd_res <- MultiCCA(lapply(X, t), ncomponents = 2, penalty = 10)
pmd_res$ws[[2]][, 1] <- -pmd_res$ws[[2]][, 1]
pmd_res$ws[[2]][, 2] <- -pmd_res$ws[[2]][, 2]
plot_simulation(mW, pmd_res)

###############################################################################
## Exact same analysis but with ordered Ws
###############################################################################
pmd_res <- MultiCCA(lapply(X, t), ncomponents = 2, type = "ordered")
p <- plot_simulation(mW, pmd_res)
ggsave("../chapter/figure/pmd/illustration_sequence.png", p[[1]], width = 5, height = 3)
ggsave("../chapter/figure/pmd/illustration_scatter.png", p[[2]], width = 5, height = 3)
