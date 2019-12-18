data {
  int<lower=0> N;
  real y[N];
  real<lower=0> sigma[N];
}

parameters {
  real mu;
  real<lower = 0.0> tau;
  real z[N];
}

transformed parameters {
  real theta[N];
  for (j in 1:N)
    theta[j] = mu + tau * z[j];
}

model {
  mu ~ normal(0, 10);
  tau ~ normal(0, 10);
  z ~ normal(0, 1);
  y ~ normal(theta, sigma);
}
