data {
  int J1;
  int J2;
  int N[J1, J2];
  int y[J1, J2];
}

parameters {
  real<lower = 0.0> sd1;
  real<lower = 0.0> sd2;
  
  real intercept;
  
  vector[J1] alpha1;
  vector[J2] alpha2;
}

model {
  sd1 ~ normal(0, 1);
  sd2 ~ normal(0, 1);
  
  intercept ~ normal(0, 1);
  
  alpha1 ~ normal(0, sd1);
  alpha2 ~ normal(0, sd2);
  
  for(i in 1:J1) {
    for(j in 1:J2) {
      y[i, j] ~ binomial_logit(N[i, j], intercept + alpha1[i] + alpha2[j]);
    }
  }
}