library(cmdstanr)
library(tidyverse)
library(posterior)

model1 = cmdstan_model("eight_schools_c.stan", quiet = FALSE)
model2 = cmdstan_model("eight_schools_nc.stan", quiet = FALSE)
model3 = cmdstan_model("eight_schools_kleppe.stan", quiet = FALSE)

data = list(N = 8,
            y = c(28,  8, -3,  7, -1,  1, 18, 12),
            sigma = c(c(15, 10, 0.16, 11, 9, 11, 10, 18)))

fit1 = model1$sample(data = data, metric = "dense_e", num_cores = 4)
fit2 = model2$sample(data = data, metric = "dense_e", num_cores = 4)
fit3 = model3$sample(data = data, metric = "dense_e", num_cores = 4)

summarise_draws(fit1$draws()) %>%
  arrange(ess_bulk) %>%
  head(2)
summarise_draws(fit2$draws()) %>%
  arrange(ess_bulk) %>%
  head(2)
summarise_draws(fit3$draws()) %>%
  arrange(ess_bulk) %>%
  head(2)


fit3$draws() %>%
  as_draws_df() %>%
  select(contains("bar")) %>%
  cor %>%
  reshape2::melt() %>%
  as_tibble() %>%
  ggplot() +
  geom_tile(aes(x = Var1, y = Var2, fill = value)) +
  coord_fixed() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab")

getEss = function(fit) {
  sdf = fit$summary()
  lp_ess = sdf %>%
    filter(variable == "lp__") %>%
    pull(ess_bulk)
  
  min_ess = sdf %>%
    filter(variable != "lp__") %>%
    arrange(ess_bulk) %>%
    head(1) %>%
    pull(ess_bulk)
  
  return(tibble(lp_ess = lp_ess, min_ess = min_ess))
}

options = list(centered = model1,
               non_centered = model2,
               kleppe = model3)
#
#seq(0.001, 1.0, length = 100)
perfdf = lapply(exp(-seq(log(0.1), log(10), length = 100)), function(scale) {
  lapply(names(options), function(option_name) {
    option = options[[option_name]]
    
    data = list(N = 8,
                y = c(28,  8, -3,  7, -1,  1, 18, 12),
                sigma = scale * c(15, 10, 16, 11,  9, 11, 10, 18))
    
    print(paste0("scale: ", scale, ", method: ", option_name))
    capture.output(fit <- option$sample(data = data,
                                        metric = "dense_e",
                                        num_cores = 4))
    
    getEss(fit) %>%
      mutate(method = option_name)
  }) %>%
    bind_rows() %>%
    mutate(scale = scale)
}) %>%
  bind_rows()

write_csv(perfdf, "eight_schools_perf_df.csv")

perfdf %>%
  gather(which_ess, ess, lp_ess, min_ess) %>%
  ggplot(aes(scale, ess)) +
  geom_line(aes(color = method)) +
  geom_point(aes(color = method)) +
  facet_grid(~ which_ess) +
  xlab("Measurement error scale (small error left, big error right)") +
  ylab("Effective sample size of either lp or the smallest non lp variable") +
  ggtitle("Varying measurement error to see parameterization performance") +
  scale_x_log10()
