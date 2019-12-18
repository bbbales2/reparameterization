data {
  int<lower=0> N;
  real y[N];
  real<lower=0> sigma[N];
}

parameters {
  real mu;
  real<lower = 0.0> tau;
  real theta[N];
}

model {
  mu ~ normal(0, 10);
  tau ~ normal(0, 10);
  theta ~ normal(mu, tau);
  y ~ normal(theta, sigma);
}
