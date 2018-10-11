#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## A simulation study using probabilistic CCA.
##
## author: sankaran.kris@gmail.com
## date: 10/04/2017

###############################################################################
## Libraries and setup
###############################################################################
library("tidyverse")
library("rstan")
source("plot.R")
set.seed(9282017)

matnorm <- function(n, p, mu = 0, sigma = 1) {
  matrix(
    rnorm(n * p, mu, sigma),
    n, p
  )
}

###############################################################################
## Simulation study to validate approach
###############################################################################
K <- 2
L1 <- 3
L2 <- 3
ps <- c(20, 50)
n <- 100

## simulate data
theta <- list(
  "Wx" = matnorm(20, L1, 0, 2),
  "Wy" = matnorm(50, L2, 0, 2),
  "Bx" = matnorm(20, K),
  "By" = matnorm(50, K),
  "xi_s" = matnorm(n, K),
  "xi_x" = matnorm(n, L1),
  "xi_y" = matnorm(n, L2)
)

x <- theta$xi_s %*% t(theta$Wx) +
  theta$xi_x %*% t(theta$Bx) +
  matnorm(n, ps[1], 0, 0.1)
y <- theta$xi_s %*% t(theta$Wy) +
  theta$xi_y %*% t(theta$By) +
  matnorm(n, ps[2], 0, 0.1)

## fit model
m <- stan_model("prob_cca.stan")
stan_data <- list(
  "K" = K,
  "L1" = L1,
  "L2" = L2,
  "p1" = ps[1],
  "p2" = ps[2],
  "n" = n,
  "tau" = 5,
  "x" = x,
  "y" = y,
  "id_x" = diag(ps[1]),
  "id_y" = diag(ps[2]),
  "id_k" = diag(K),
  "zeros_k" = rep(0, K)
)

fit <- vb(m, stan_data, iter = 5000)
theta_samples <- rstan::extract(fit)
rm(fit)
theta_hat <- parameter_means(theta_samples)
plot(theta_hat$xi_s[, 1], theta$xi_s[, 2], asp = 1)
