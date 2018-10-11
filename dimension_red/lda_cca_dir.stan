/*
 * Fit a LDA + CCA model
 *
 * author: sankaran.kris@gmail.com
 * date: 10/03/2017
 */

data {
  int<lower=1> n;
  int<lower=1> p1;
  int<lower=1> p2;
  int<lower=1> K;
  int<lower=1> L1;
  int<lower=1> L2;
  real<lower=0> sigma;
  vector<lower=0>[K] alpha_k;
  vector<lower=0>[L1] alpha_l1;
  vector<lower=0>[L2] alpha_l2;
  vector<lower=0>[p1] gamma;

  int<lower=0> x[n, p1];
  matrix[n, p2] y;
  matrix[p2, p2] id_y;
}

parameters {
  simplex[K] xi_s[n];
  simplex[L1] xi_x[n];
  simplex[L2] xi_y[n];
  matrix[L1, p1] Wx;
  matrix[L2, p2] Wy;
  matrix[K, p1] Bx;
  matrix[K, p2] By;
}

transformed parameters {
  matrix[p2, p2] cov_y;
  cov_y = sigma ^ 2 * id_y;
}

model {
  // prior
  for (i in 1:n) {
    xi_s[i] ~ dirichlet(alpha_k);
    xi_x[i] ~ dirichlet(alpha_l1);
    xi_y[i] ~ dirichlet(alpha_l2);
  }

  for (k in 1:K) {
    Bx[k]' ~ dirichlet(gamma);
  }

  for (l in 1:L1) {
    Wx[l]' ~ dirichlet(gamma);
  }

  // likelihood
  for (i in 1:n) {
    x[i] ~ multinomial(Bx' * xi_s[i] + Wx' * xi_x[i]);
    y[i] ~ multi_normal(By' * xi_s[i] + Wy' * xi_y[i], cov_y);
  }
}
