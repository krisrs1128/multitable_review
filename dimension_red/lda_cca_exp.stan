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
  real<lower=0> a0;
  real<lower=0> b0;

  int<lower=0> x[n, p1];
  matrix[n, p2] y;
  matrix[p2, p2] id_y;
  matrix[K, K] id_k;
  matrix[L1, L1] id_l1;
  matrix[L2, L2] id_l2;
  vector<lower=0, upper=0>[K] zeros_k;
  vector<lower=0, upper=0>[L1] zeros_l1;
  vector<lower=0, upper=0>[L2] zeros_l2;
}


parameters {
  real<lower=0> tau_sq[3];
  matrix[n, K] xi_s;
  matrix[n, L1] xi_x;
  matrix[n, L2] xi_y;
  matrix<lower=0>[p1, L1] Wx;
  matrix<lower=0>[p2, L2] Wy;
  matrix<lower=0>[p1, K] Bx;
  matrix<lower=0>[p2, K] By;
}

transformed parameters {
  matrix[p2, p2] cov_y;
  real<lower=0> tau[3];

  cov_y = sigma ^ 2 * id_y;
  tau = sqrt(tau_sq);
}

model {
  // prior
  for (v in 1:3) {
    tau_sq[v] ~ inv_gamma(a0, b0);
  }

  for (i in 1:n) {
    xi_s[i] ~ multi_normal(zeros_k, tau[1] * id_k);
    xi_x[i] ~ multi_normal(zeros_l1, tau[2] * id_l1);
    xi_y[i] ~ multi_normal(zeros_l2, tau[3] * id_l2);
  }

  for (j in 1:p1) {
    for (k in 1:K) {
      Bx[j, k] ~ double_exponential(0, 1);
    }
    for (l in 1:L2) {
      Wx[j, l] ~ double_exponential(0, 1);
    }
  }

  for (j in 1:p2) {
    for (k in 1:K) {
      By[j, k] ~ double_exponential(0, 1);
    }
    for (l in 1:L2) {
      Wy[j, l] ~ double_exponential(0, 1);
    }
  }

  // likelihood
  for (i in 1:n) {
    x[i] ~ multinomial(softmax(Bx * xi_s[i]' + Wx * xi_x[i]'));
    y[i] ~ multi_normal(By * xi_s[i]' + Wy * xi_y[i]', cov_y);
  }
}
