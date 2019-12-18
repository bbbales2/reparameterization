data {
  int J1;
  int J2;
  int N[J1, J2];
  int y[J1, J2];
  
  int Fsd;
  real sd1f[Fsd];
  real sd2f[Fsd];
  
  int Falpha;
  real interceptf[Falpha];
  vector[J1] alpha1f[Falpha];
  vector[J2] alpha2f[Falpha];
}

parameters {
  real<lower = 0.0> sd1[1 - Fsd];
  real<lower = 0.0> sd2[1 - Fsd];
  
  real intercept[1 - Falpha];
  
  vector[J1] alpha1[1 - Falpha];
  vector[J2] alpha2[1 - Falpha];
}

model {
  real sd1a = (Fsd) ? sd1f[1] : sd1[1];
  real sd2a = (Fsd) ? sd2f[1] : sd2[1];
  
  real intercepta = (Falpha) ? interceptf[1] : intercept[1];
  
  vector[J1] alpha1a = (Falpha) ? alpha1f[1] : alpha1[1];
  vector[J2] alpha2a = (Falpha) ? alpha2f[1] : alpha2[1];
  
  sd1a ~ normal(0, 1);
  sd2a ~ normal(0, 1);
  
  intercepta ~ normal(0, 1);
  
  alpha1a ~ normal(0, sd1a);
  alpha2a ~ normal(0, sd2a);
  
  if(Falpha == 0) {
    for(i in 1:J1) {
      for(j in 1:J2) {
        y[i, j] ~ binomial_logit(N[i, j], intercepta + alpha1a[i] + alpha2a[j]);
      }
    }
  }
}