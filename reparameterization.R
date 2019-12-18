library(cmdstanr)
library(rstan)
library(tidyverse)
library(posterior)

set_cmdstan_path("/home/bbales2/cmdstan-latest")

J1 = 4
J2 = 10

sd1 = 0.3
sd2 = 0.15

alpha1 = rnorm(J1, 0, sd1)
alpha2 = rnorm(J2, 0, sd2)

#rse = function(n, p, sel, seh) {
#  if(p <= 0.5) {
#    return(exp(runif(n, min = log(sel), max = log(sel) + p * 2 * (log(seh) - log(sel)))))
    #return(c(min = sel, max = sel + p * 2 * (seh - sel)))
#  } else {
#    return(exp(runif(n, min = log(sel) + ((p - 0.5) * 2) * (log(seh) - log(sel)), max = log(seh))))
    #return(c(min = sel + ((p - 0.5) * 2) * (seh - sel), max = seh))
#  }
#}

rse = function(n, x, sel, seh) {
  logx = rep(0, n)
  
  logx = rbeta(n, 9 * x^2 + 1, 9 * abs((x - 1.0)^2) + 1)

  exp(log(sel) + logx * (log(seh) - log(sel)))
}

intercept = 0.25

inv_logit = function(x) 1 / (1 + exp(-x))

generate_data = function(x, sel, seh) {
  se = matrix(rse(J1 * J2, x, sel, seh), nrow = J1)
  
  p = array(0, dim = c(J1, J2))
  Nj = array(0, dim = c(J1, J2))
  y = array(0, dim = c(J1, J2))
  for(i in 1:J1) {
    for(j in 1:J2) {
      p[i, j] = inv_logit(intercept + alpha1[i] + alpha2[j])
      Nj[i, j] = ceiling(p[i, j] * (1 - p[i, j]) / se[i, j]^2)
      y[i, j] = rbinom(1, Nj[i, j], p[i, j])
    }
  }
  
  return(list(Nj = Nj, y = y, p = round(p, 2), se = round(se, 5)))
}

data_est = list(J1 = J1, J2 = J2, N = Nj, y = y,
                Fsd = 0, sd1f = array(0, dim = c(0)), sd2f = array(0, dim = c(0)),
                Falpha = 0, interceptf = array(0, dim = c(0)), alpha1f = array(0, dim = c(0, J1)), alpha2f = array(0, dim = c(0, J2)))
data_fix_sds = function(sd1, sd2) {
  list(J1 = J1, J2 = J2, N = Nj, y = y,
       Fsd = 1, sd1f = array(sd1, dim = c(1)), sd2f = array(sd2, dim = c(1)),
       Falpha = 0, interceptf = array(0, dim = c(0)), alpha1f = array(0, dim = c(0, J1)), alpha2f = array(0, dim = c(0, J2)))
}
data_fix_alpha = function(intercept, alpha1, alpha2) {
  list(J1 = J1, J2 = J2, N = Nj, y = y,
       E = 0, sd1f = array(sd1, dim = c(1)), sd2f = array(sd2, dim = c(1)),
       Falpha = 1, interceptf = array(intercept, dim = c(1)), alpha1f = array(alpha1, dim = c(1, J1)), alpha2f = array(alpha2, dim = c(1, J2)))
}

models = stan_model("regression.stan")
o = optimizing(models, data = data_fix_sds(1.0, 1.0), hessian = TRUE)

modelc = cmdstan_model("regression_c.stan", quiet = FALSE)
modeln = cmdstan_model("regression_nc.stan", quiet = FALSE)
modelk = cmdstan_model("regression_kleppe.stan", quiet = FALSE)

fitc = modelc$sample(data = data_est, metric = "dense_e", num_cores = 4)
fitn = modeln$sample(data = list(J1 = J1, J2 = J2, N = Nj, y = y), metric = "dense_e", num_cores = 4)
fitk = modelk$sample(data = list(J1 = J1, J2 = J2, N = Nj, y = y,
                                 intercept_mu = o$par[1],
                                 alpha1_mu = o$par[2:(1 + J1)],
                                 alpha2_mu = o$par[(2 + J1):(1 + J1 + J2)],
                                 L = t(chol(-o$hessian))), metric = "dense_e", num_cores = 4)

fit$summary() %>% arrange(ess_bulk) %>% head(2)
fitn$summary() %>% arrange(ess_bulk) %>% head(2)
fitk$summary() %>% arrange(ess_bulk) %>% head(2)

fitk$draws() %>%
  posterior::as_draws_df() %>%
  select(starts_with("mu"), .chain, .iteration) %>%
  posterior::summarise_draws()

fit$draws() %>%
  as_draws_df() %>%
  select(starts_with("intercept"), starts_with("alpha1"), starts_with("alpha2"), .chain, .iteration) %>%
  summarise_draws()

methods = list(centered = modelc,
               noncentered = modeln,
               transformed = modelk)

perfdf = lapply(seq(0.0, 1.0, length = 100), function(scale) {
  data = generate_data(scale, 0.01, 10.0)
  data_fix_sds = function(sd1, sd2) {
    list(J1 = J1, J2 = J2, N = data$Nj, y = data$y,
         Fsd = 1, sd1f = array(sd1, dim = c(1)), sd2f = array(sd2, dim = c(1)),
         Falpha = 0, interceptf = array(0, dim = c(0)), alpha1f = array(0, dim = c(0, J1)), alpha2f = array(0, dim = c(0, J2)))
  }
  o = optimizing(models, data = data_fix_sds(1.0, 1.0), hessian = TRUE)
  lapply(names(methods), function(method_name) {
    print(paste0("scale: ", scale, ", method: ", method_name))
    
    capture.output(fit <- methods[[method_name]]$sample(data = list(J1 = J1, J2 = J2, N = data$Nj, y = data$y,
                                                    intercept_mu = o$par[1],
                                                    alpha1_mu = o$par[2:(1 + J1)],
                                                    alpha2_mu = o$par[(2 + J1):(1 + J1 + J2)],
                                                    L = t(chol(-o$hessian))), metric = "dense_e", num_cores = 4))
  
    min_ess = fit$draws() %>%
      as_draws_df %>%
      select(starts_with("sd"), starts_with("intercept"), starts_with("alpha"), -contains("bar"), -contains("z"), .chain, .iteration) %>%
      summary() %>%
      arrange(ess_bulk) %>%
      head(1) %>%
      pull(ess_bulk)
    
    min_lp_ess = fit$summary() %>%
      filter(variable == "lp__") %>%
      pull(ess_bulk)
    
    return(tibble(method = method_name, min_ess = min_ess, min_lp_ess = min_lp_ess))
  }) %>% bind_rows() %>%
    mutate(scale = scale)
}) %>% bind_rows()

perfdf %>%
  gather(which_ess, ess, min_lp_ess, min_ess) %>%
  ggplot(aes(scale, ess)) +
  geom_line(aes(color = method)) +
  geom_point(aes(color = method)) +
  facet_grid(~ which_ess) +
  xlab("Measurement error scale (small error left, big error right)") +
  ylab("Smallest effective sample size in model (excluding lp)")
#+
 # scale_x_log10()
