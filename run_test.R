i <- 13
dat <- list()
dat$cases <- deaths_community[age_grp == agelabs[i], confirm][40:125]
dat$date <- deaths_community[age_grp == agelabs[i], date][40:125]
dat$t <- length(dat$cases)
dat$time <- 1:dat$t
dat$M <- ceiling(0.3 * dat$t)
dat$L <- dat$t * 2
dat$delays <- 2
dat$delay_mean_sd <- c(delays[[i]]$mean_sd, incubation_period$mean_sd)
dat$delay_mean_mean <- c(delays[[i]]$mean, incubation_period$mean)
dat$delay_sd_mean <- c(delays[[i]]$sd, incubation_period$sd)
dat$delay_sd_sd <- c(delays[[i]]$sd_sd, incubation_period$sd_sd)
dat$max_delay <- c(delays[[i]]$max, incubation_period$max)
dat$gt_mean_mean <- generation_time$mean
dat$gt_mean_sd <- generation_time$mean_sd
dat$gt_sd_mean <- generation_time$sd
dat$gt_sd_sd <- generation_time$sd_sd
dat$max_gt <- generation_time$max
dat$lengthscale_alpha <- 4.5
dat$lengthscale_beta <- 21.5
dat$alpha_sd <- 0.1

mod <- rstan::stan_model(here::here("test.stan"))
seedx <- as.integer(paste( sample( 1:9, 8, replace=TRUE ), collapse="" ))
fit <- rstan::sampling(mod, 
                       data = dat, 
                       chains = 1, 
                       seed = seedx,
                       control = list(max_treedepth = 12, adapt_delta = 0.9))
res <- rstan::extract(fit)
pairs(fit, pars = c("rho", "alpha", "phi", "lp__", "log_initial_inf"))
pairs(fit, pars = c("log_initial_inf", paste0("log_R[", 5:10, "]"), "lp__"))
pairs(fit, pars = c("log_initial_inf", "lp__"))
pairs(fit, pars = c("gt_mean", "gt_sd", "delay_mean", "delay_sd", "lp__"))
pairs(fit, pars = c(paste0("beta[", 21:25, "]"), "lp__"))

draws <- as.array(fit)
np <- nuts_params(fit)

color_scheme_set("darkgray")
div_style <- parcoord_style_np(div_color = "green", div_size = 0.5, div_alpha = 1)
mcmc_parcoord(draws, np = np, pars = c(paste0("log_R[", 1:86, "]")), np_style = div_style)
mcmc_parcoord(draws, np = np, pars = c(paste0("R[", 1:86, "]")), np_style = div_style)

R_tab <- melt(data.table(res$R))[, t := rep(1:ncol(res$R), rep(1000, ncol(res$R)))
                                 ][, iter := rep(1:1000, ncol(res$R))
                                   ][, variable := NULL
                                     ][, date := rep(dat$date, rep(1000, length(dat$date)))]

inf_tab <- melt(data.table(res$infections))[, t := rep(1:ncol(res$infections), rep(1000, ncol(res$infections)))
][, iter := rep(1:1000, ncol(res$infections))
][, variable := NULL][, date := rep(dat$date, rep(1000, length(dat$date)))]

rep_tab <- melt(data.table(res$reports))[, t := rep(1:ncol(res$reports), rep(1000, ncol(res$reports)))
][, iter := rep(1:1000, ncol(res$reports))
][, variable := NULL][, date := rep(dat$date, rep(1000, length(dat$date)))]


p1 <- R_tab[, .(med = median(value),
          top = quantile(value, prob = 0.975),
          bottom = quantile(value, prob = 0.025)), by = "date"] %>%
  ggplot(aes(x = as.Date(date), y = med)) +
  geom_ribbon(aes(ymax = top, ymin = bottom), alpha = 0.5) + 
  geom_line() +
  labs(y = "Rt")

p2 <- rep_tab[, .(med = median(value),
          top = quantile(value, prob = 0.975),
          bottom = quantile(value, prob = 0.025)), by = "date"] %>%
  ggplot(aes(x = as.Date(date), y = med)) +
  geom_ribbon(aes(ymax = top, ymin = bottom), alpha = 0.5) + 
  geom_line() +
  geom_line(data = data.table(t = dat$date, med = dat$cases), inherit.aes = FALSE, aes(x = t, y= med), col = "red") +
  labs(y = "Reported Deaths")


p3 <- inf_tab[, .(med = median(value),
                top = quantile(value, prob = 0.975),
                bottom = quantile(value, prob = 0.025)), by = "date"] %>%
  ggplot(aes(x = as.Date(date), y = med)) +
  geom_ribbon(aes(ymax = top, ymin = bottom), alpha = 0.5) + 
  geom_line() +
  labs(y = "Infections") 

p1 / p2 / p3  &
  geom_vline(xintercept = as.Date("2020-03-16"), lty = 2) & labs(x = "Date")


