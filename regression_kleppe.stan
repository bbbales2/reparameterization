data {
  int J1;
  int J2;
  int N[J1, J2];
  int y[J1, J2];

  real intercept_mu;
  vector[J1] alpha1_mu;
  vector[J2] alpha2_mu;
  
  matrix[1 + J1 + J2, 1 + J1 + J2] L;
}

parameters {
  real<lower = 0.0> sd1;
  real<lower = 0.0> sd2;
  
  real intercept_bar;
  
  vector[J1] alpha1_bar;
  vector[J2] alpha2_bar;
}

transformed parameters {
  real intercept;
  
  vector[J1] alpha1;
  vector[J2] alpha2;
  
  vector[1 + J1 + J2] mu;

  real adj_jac = 0.0;

  {
    vector[1 + J1 + J2] beta_bar = append_row(intercept_bar, append_row(alpha1_bar, alpha2_bar));
    vector[1 + J1 + J2] beta;
    vector[1 + J1 + J2] mu_bar = append_row(intercept_mu, append_row(alpha1_mu, alpha2_mu));
    
    matrix[1 + J1 + J2, 1 + J1 + J2] KL = L;

    KL[1, 1] = KL[1, 1] + 1;
    
    for(i in 1:J1) {
      KL[1 + i, 1 + i] = KL[1 + i, 1 + i] + 1 / sd1;
    }

    for(i in 1:J2) {
      KL[1 + J1 + i, 1 + J1 + i] = KL[1 + J1 + i, 1 + J1 + i] + 1 / sd2;
    }

    for(i in 1:(1 + J1 + J2)) {
      adj_jac = adj_jac - log(KL[i, i]);
    }
    
    {
      vector[1 + J1 + J2] zmu = L * (L' * mu_bar);
      vector[1 + J1 + J2] pmu = mdivide_left_tri_low(KL, zmu);
      mu = mdivide_right_tri_low(pmu', KL)';
    }
      
    beta = mdivide_right_tri_low(beta_bar', KL)';
    
    intercept = mu[1] + beta[1];
    alpha1 = mu[2 : (1 + J1)] + beta[2 : (1 + J1)];
    alpha2 = mu[(2 + J1) : (1 + J1 + J2)] + beta[(2 + J1) : (1 + J1 + J2)];
  }
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
  
  target += adj_jac;
}