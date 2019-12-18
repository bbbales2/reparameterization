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
  real jac_adj = 0.0;
  
  for(n in 1:N) {
    real L = sqrt(1 / sigma[n]^2 + 1 / tau^2);
    theta[n] = (mu / tau^2 + y[n] / sigma[n]^2) / L^2 + z[n] / L;
    jac_adj = jac_adj - log(L);
  }
}

model {
  mu ~ normal(0, 10);
  tau ~ normal(0, 10);
  theta ~ normal(mu, tau);
  y ~ normal(theta, sigma);
  
  target += jac_adj;
}
