#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## An experiment using Latent Dirichlet Allocation + Canonical Correlation
## Analysis simultaneously across count and gaussian tables.
##
## author: sankaran.kris@gmail.com
## date: 10/03/2017

###############################################################################
## Libraries and setup
###############################################################################
library("rstan")
library("tidyverse")
library("reshape2")
set.seed(1032017)

#' i.i.d. Gaussian Matrix
matnorm <- function(n, p, mu = 0, sigma = 1) {
  matrix(
    rnorm(n * p, mu, sigma),
    n, p
  )
}

#' Simulate parameters used in LDA + CCA model
simulate_parameters <- function(n = 100, N = 1000, p1 = 200, p2 = 30, K = 2,
                                L1 = 3, L2 = 3) {
  list(
    "n" = n,
    "p1" = p1,
    "p2" = p2,
    "K" = as.integer(K),
    "L1" = as.integer(L1),
    "L2" = as.integer(L2),
    "N" = 1000,
    "xi_s" = matnorm(n, K),
    "xi_x" = matnorm(n, L1),
    "xi_y" = matnorm(n, L2),
    "Bx" = matnorm(p1, K),
    "By" = matnorm(p2, K),
    "Wx" = matnorm(p1, L2),
    "Wy" = matnorm(p2, L1),
    "sigma" = 0.2
  )
}

#' Softmax function
softmax <- function(x) {
  exp(x) / sum(exp(x))
}

#' Simulate from LDA + CCA model
simulate_data <- function(theta) {
  mu_x <- theta$xi_s %*% t(theta$Bx) + theta$xi_x %*% t(theta$Wx)
  mu_y <- theta$xi_s %*% t(theta$By) + theta$xi_y %*% t(theta$Wy)

  X <- matrix(nrow = theta$n, ncol = theta$p1)
  Y <- matrix(nrow = theta$n, ncol = theta$p2)

  for (i in seq_len(theta$n)) {
    X[i, ] <- rmultinom(1, theta$N, softmax(mu_x[i, ]))
  }

  for (i in seq_len(theta$n)) {
    for (j in seq_len(theta$p2)) {
      Y[i, j] <- rnorm(1, mu_y[i, j], theta$sigma)
    }
  }

  list("X" = X, "Y" = Y)
}

###############################################################################
## simulate toy data
###############################################################################
theta <- simulate_parameters()
sim <- simulate_data(theta)
plot(sim$X[, 121], sim$Y[, 11])
cormat <- cor(sim$X, sim$Y)
heatmap(cormat)

###############################################################################
## use Rstan model
###############################################################################

stan_data <- list(
  "n" = theta$n,
  "p1" = theta$p1,
  "p2" = theta$p2,
  "K" = theta$K,
  "L1" = theta$L1,
  "L2" = theta$L1,
  "sigma" = theta$sigma,
  "tau" = 5,
  "x" = sim$X,
  "y" = sim$Y,
  "id_y" = diag(theta$p2),
  "id_k" = diag(theta$K),
  "id_l1" = diag(theta$L1),
  "id_l2" = diag(theta$L2),
  "zeros_k" = rep(0, theta$K),
  "zeros_l1" = rep(0, theta$L1),
  "zeros_l2" = rep(0, theta$L2)
)

m <- stan_model("lda_cca.stan")
stan_fit <- vb(m, data = stan_data)


###############################################################################
## Inspect results
###############################################################################
samples <- rstan::extract(stan_fit)
xi_s_hat <- apply(samples$xi_s, c(2, 3), mean)
Bx_hat <- apply(samples$Bx, c(2, 3), mean)

## seems to work!!
plot(xi_s_hat[, 1], theta$xi_s[, 1])
plot(xi_s_hat[, 2], theta$xi_s[, 2])
plot(Bx_hat[, 2], theta$Bx[, 2])
