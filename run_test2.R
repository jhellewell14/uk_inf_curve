example_confirmed <- EpiNow2::example_confirmed

reporting_delay <- bootstrapped_dist_fit(rlnorm(100, log(4), 1), max_value = 30)
generation_time <- get_generation_time(disease = "SARS-CoV-2", source = "ganyani")
incubation_period <- get_incubation_period(disease = "SARS-CoV-2", source = "lauer")

dat <- list()
# dat$cases <- deaths_community[age_grp == agelabs[i], confirm][40:260]
# dat$date <- deaths_community[age_grp == agelabs[i], date][40:260]
dat$cases <- c(rep(0, 7), example_confirmed[, confirm])
# dat$cases <- example_confirmed[1:80, confirm]
dat$date <- c(seq(from = min(example_confirmed$date) - 7, to = min(example_confirmed$date) - 1, by = "day"), 
	example_confirmed[, date])
# dat$date <- example_confirmed[1:80, date]
dat$day_of_week <- lubridate::wday(dat$date, week_start = 1)
dat$t <- length(dat$cases)
# dat$n <- 10
# dat$date_n <- tab$date[(60 - dat$n): 292]
dat$time <- 1:(dat$t)
dat$M <- ceiling(0.3 * dat$t)
dat$L <- dat$t * 2
dat$delays <- 2
dat$delay_mean_sd <- c(reporting_delay$mean_sd, incubation_period$mean_sd)
dat$delay_mean_mean <- c(reporting_delay$mean, incubation_period$mean)
dat$delay_sd_mean <- c(reporting_delay$sd, incubation_period$sd)
dat$delay_sd_sd <- c(reporting_delay$sd_sd, incubation_period$sd_sd)
dat$max_delay <- c(reporting_delay$max, incubation_period$max)
dat$gt_mean_mean <- generation_time$mean
dat$gt_mean_sd <- generation_time$mean_sd
dat$gt_sd_mean <- generation_time$sd
dat$gt_sd_sd <- generation_time$sd_sd
dat$max_gt <- generation_time$max
dat$lengthscale_alpha <- 4.5
dat$lengthscale_beta <- 21.5
dat$alpha_sd <- 2


mod <- rstan::stan_model("~/repos/uk_inf_curve/test.stan")
fit <- rstan::sampling(mod,
                       data = dat,
                       chains = 1,
                       seed = 54321,
                       control = list(max_treedepth = 12, adapt_delta = 0.8))
res <- rstan::extract(fit)


library(bayesplot)

draws <- as.array(fit, pars = paste0("log_initial_inf[",8:137,"]"))
np <- nuts_params(fit)
p1 <- mcmc_parcoord(draws, alpha = 0.05) + ggtitle("Log infections (GP)")


draws <- as.array(fit, pars = paste0("infections[",8:137,"]"))
np <- nuts_params(fit)
p3 <- mcmc_parcoord(draws, alpha = 0.05) + ggtitle("Infections")

draws <- as.array(fit, pars = paste0("R[",8:137,"]"))
np <- nuts_params(fit)
p2 <- mcmc_parcoord(draws, alpha = 0.05) + ggtitle("R")

out <- p1 / p2 / p3 

reported_cases <- example_confirmed

estimates <- epinow(reported_cases = reported_cases, 
                    generation_time = generation_time,
                    delays = list(incubation_period, reporting_delay))

outy <- estimates$plot

outy[[5]]
ggsave(out, filename = "out3.png")

