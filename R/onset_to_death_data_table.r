library(data.table)
library(magrittr)
# Set number of threads for data.table
setDTthreads(parallel::detectCores())

# Read in data
data <- data.table::fread("~/Downloads/CCPUKSARI_DATA_2020-08-04_0947.csv", na.strings = "")

# Select columns + fix read in issue where entries are "" instead of NA
df <- data[,.(cestdat, dsstdtc, dsterm, subjid, age = age_estimateyears)]


df[, c("onset_date_missing", "outcome_date_missing", "dead") :=
      list(all(is.na(cestdat)), all(is.na(dsstdtc)), any(dsterm == 4, na.rm = TRUE)), 
   by = "subjid"]

df <- df[!onset_date_missing & !outcome_date_missing & dead
][, .(onset_date = unique(cestdat[!is.na(cestdat)]),
      dead = unique(dead),
      age = unique(age[!is.na(age)]),
      outcome_date = unique(dsstdtc[!is.na(dsstdtc)])), by = "subjid"
][, delay := as.integer(as.Date(outcome_date) - as.Date(onset_date))
][delay >= 0 & delay <= 60 & as.Date(onset_date) > "2020-01-01"][
   , delay_sampled := runif(.N, delay, delay + 1)]

df[, age_grp := cut(age, breaks = agebreaks, labels = agelabs, right = FALSE)]

# Fit a gamma distribution
nbfit30 <- fitdistrplus::fitdist(df[age_grp == "30-39", delay_sampled], distr = "gamma")
nbfit40 <- fitdistrplus::fitdist(df[age_grp == "40-49", delay_sampled], distr = "gamma")
nbfit50 <- fitdistrplus::fitdist(df[age_grp == "50-59", delay_sampled], distr = "gamma")
nbfit60 <- fitdistrplus::fitdist(df[age_grp == "60-69", delay_sampled], distr = "gamma")
nbfit70 <- fitdistrplus::fitdist(df[age_grp == "70-79", delay_sampled], distr = "gamma")
nbfit80 <- fitdistrplus::fitdist(df[age_grp == "80-89", delay_sampled], distr = "gamma")
nbfit90 <- fitdistrplus::fitdist(df[age_grp == "90-99", delay_sampled], distr = "gamma")

y <- c(dgamma(seq(0, 60, 0.1), shape = nbfit30$estimate[1], rate = nbfit30$estimate[2]),
       dgamma(seq(0, 60, 0.1), shape = nbfit40$estimate[1], rate = nbfit40$estimate[2]),
       dgamma(seq(0, 60, 0.1), shape = nbfit50$estimate[1], rate = nbfit50$estimate[2]),
       dgamma(seq(0, 60, 0.1), shape = nbfit60$estimate[1], rate = nbfit60$estimate[2]),
       dgamma(seq(0, 60, 0.1), shape = nbfit70$estimate[1], rate = nbfit70$estimate[2]),
       dgamma(seq(0, 60, 0.1), shape = nbfit80$estimate[1], rate = nbfit80$estimate[2]),
       dgamma(seq(0, 60, 0.1), shape = nbfit90$estimate[1], rate = nbfit90$estimate[2]))

hi <- data.frame(y,
                 x = rep(seq(0, 60, 0.1), 7),
                 agegrp = rep(agelabs[4:10], rep(601, 7)))


hi %>% 
   ggplot2::ggplot(ggplot2::aes(x = x, y = y, col = as.factor(agegrp))) +
   ggplot2::geom_line() + 
   cowplot::theme_cowplot() +
   ggplot2::scale_color_discrete(name = "Age group") +
   ggplot2::labs(x = "Days since onset", y = "Probability density")


# Death distribution for care homes in linelist from covid19_automation

path_to_factory <- "~/repos/covid19_automation"

file_path <- file.path(path_to_factory, "data", "rds", "deaths_eng_latest.rds")
key <- cyphr::data_key(file.path(path_to_factory, "data"))
x <- cyphr::decrypt(readRDS(file_path), key)

lldf <- data.table::as.data.table(x)
lldf[, care_home_death := fifelse(residence_type == "care_nursing_home" |
                                     place_of_death == "care_home",
                                  "Care home",
                                  "Other")]

### Age proportions

agebreaks <- c(0, 10, 20, 30, 40, 50, 60, 70 ,80, 90, 100)
agelabs <- c("0-9", "10-19", 
             "20-29", "30-39", 
             "40-49", "50-59", 
             "60-69", "70-79",
             "80-89","90-99")

age_props <- lldf[, age_grp := cut(age, breaks = agebreaks, labels = agelabs, right = FALSE)
][, .N, age_grp][
   !is.na(age_grp)][
      , prop := N/sum(N)]


age_props[order(N)] %>%
   ggplot2::ggplot(ggplot2::aes(x = age_grp, y = prop)) +
   ggplot2::geom_bar(stat = "identity") +
   ggplot2::scale_y_continuous(breaks = seq(0, 0.4, 0.1), labels = seq(0, 40, 10)) +
   ggplot2::labs(y = "Proportion of deaths (%)", x = "Age group") +
   cowplot::theme_cowplot()

####

lldf <- lldf[!is.na(date_onset) & !is.na(date_death) & care_home_death == "Care home", delay := as.numeric(date_death - date_onset)
][delay > 0 & delay < 60 & date_onset > "2020-01-01", "delay"
]

lldf[, delay_sampled := runif(.N, delay, delay + 1)]

nbfit2 <- fitdistrplus::fitdist(lldf$delay_sampled, distr = "gamma")


# Plot fitted distribution
plot(seq(0, 60, 0.1), dgamma(seq(0, 60, 0.1), shape = nbfit$estimate[1], rate = nbfit$estimate[2]), 
     type = "l",
     ylab = "Density",
     xlab = "Days",
     main = "Onset to death delay",
     ylim = c(0, 0.1))
lines(seq(0, 60, 0.1), dgamma(seq(0, 60, 0.1), shape = nbfit2$estimate[1], rate = nbfit2$estimate[2]),
      col = "red")
text(15, 0.06, "Care homes", col = "red")
text(25, 0.03, "Hospital")


#### Deaths time series
deaths_ts <- data.table::as.data.table(x)
deaths_ts[, care_home_death := fifelse(residence_type == "care_nursing_home" |
                                          place_of_death == "care_home",
                                       "Care home",
                                       "Other")]


reported_cases_community  <- deaths_ts[ons == "reported_by_ons" & care_home_death == "Other", "date_death"][, .N, by = "date_death"][, .(confirm = N, date = date_death)][order(date)][.(seq.Date(from = min(date), to = max(date), by = "day")), 
                                                                                                                                                                                       on = .(date),roll = 0][is.na(confirm), confirm := 0]

reported_cases_community %>% 
   ggplot2::ggplot(ggplot2::aes(x = date, y = confirm)) +
   ggplot2::geom_line()

### EpiNow 2 fit

generation_time <- list(mean = EpiNow2::covid_generation_times[1, ]$mean,
                        mean_sd = EpiNow2::covid_generation_times[1, ]$mean_sd,
                        sd = EpiNow2::covid_generation_times[1, ]$sd,
                        sd_sd = EpiNow2::covid_generation_times[1, ]$sd_sd,
                        max = 30)

incubation_period <- list(mean = EpiNow2::covid_incubation_period[1, ]$mean,
                          mean_sd = EpiNow2::covid_incubation_period[1, ]$mean_sd,
                          sd = EpiNow2::covid_incubation_period[1, ]$sd,
                          sd_sd = EpiNow2::covid_incubation_period[1, ]$sd_sd,
                          max = 30)

reporting_delay <- EpiNow2::bootstrapped_dist_fit(values = df$delay, verbose = TRUE)
## Set max allowed delay to 30 days to truncate computation
reporting_delay$max <- 60


estimates <- EpiNow2::estimate_infections(reported_cases = reported_cases_community, 
                                          generation_time = generation_time,
                                          estimate_rt = FALSE, fixed = FALSE,
                                          delays = list(incubation_period, reporting_delay),
                                          horizon = 7, samples = 4000, warmup = 500, 
                                          cores = 4, chains = 4, verbose = TRUE, 
                                          adapt_delta = 0.95)

### Plot

p1 <- estimates$summarised[variable == "infections" & type == "estimate",] %>%
   ggplot2::ggplot(ggplot2::aes(x = date, y= median)) +
   ggplot2::geom_line(col = "dodgerblue") +
   ggplot2::geom_ribbon(data = reported_cases_community, ggplot2::aes(x = date, ymax = confirm, ymin = 0), 
                        inherit.aes = FALSE, lty = 2, fill = "red4", alpha = 0.8) +
   ggplot2::geom_ribbon(ggplot2::aes(ymin = bottom, ymax = top), alpha = 0.25, fill = "dodgerblue") +
   ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.25, fill = "dodgerblue") +
   ggplot2::geom_vline(xintercept = as.Date("2020-03-23"), lty = 2) +
   cowplot::theme_cowplot() +
   ggplot2::labs(y = "Daily incidence", x = "Date")

### Proportion of hospital infections in non-carehome data

df2 <- data[,.(cestdat, dsstdtc, dsterm, subjid, hostdat)]

df2 <- df2[, c("onset_date_missing", "outcome_date_missing", "hosp_date_missing", "dead") :=
              list(all(is.na(cestdat)), all(is.na(dsstdtc)), all(is.na(hostdat)), any(dsterm == 4, na.rm = TRUE)), 
           by = "subjid"][!onset_date_missing & !outcome_date_missing & !hosp_date_missing & dead
           ][, .(onset_date = unique(cestdat[!is.na(cestdat)]),
                 dead = unique(dead),
                 outcome_date = unique(dsstdtc[!is.na(dsstdtc)]),
                 hosp_date = unique(hostdat[!is.na(hostdat)])), by = "subjid"
           ]


df2 <- df2[order(outcome_date), .(hosp = sum((as.Date(onset_date) >= as.Date(hosp_date) + 5), na.rm = TRUE), 
                                  N = .N), by = "outcome_date"][, prop := hosp/N, ] 

df2 %>%
   ggplot2::ggplot(ggplot2::aes(x = outcome_date)) +
   ggplot2::geom_ribbon(ggplot2::aes( ymin = 0, ymax = N), fill = "black") +
   ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = hosp), fill = "yellow") +
   cowplot::theme_cowplot() +
   ggplot2::labs(x = "Date", y = "Deaths")

df2[, ind := 1:.N] %>%
   ggplot2::ggplot(ggplot2::aes(x = ind, y = prop, col = N)) +
   ggplot2::geom_point() +
   cowplot::theme_cowplot() +
   ggplot2::labs(x = "Time", y = "Proportion of deaths from hospital-acquired infections (%)") +
   ggplot2::scale_color_continuous(name = "Total deaths") +
   ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.1), labels = seq(0, 100, 10)) +
   ggplot2::geom_hline(yintercept = median(df2$prop), lty = 2)

### Reported deaths in care homes
reported_cases_carehome <- reported_cases_community  <- deaths_ts[ons == "reported_by_ons" & care_home_death == "Care home", "date_death"
][, .N, by = "date_death"
][, .(confirm = N, date = date_death)
][order(date)
][.(seq.Date(from = min(date), to = max(date), by = "day")), on = .(date), roll = 0
][is.na(confirm), confirm := 0]

reporting_delay_ch <- EpiNow2::bootstrapped_dist_fit(values = lldf$delay, verbose = TRUE)
## Set max allowed delay to 30 days to truncate computation
reporting_delay_ch$max <- 60

estimates_ch <- EpiNow2::estimate_infections(reported_cases = reported_cases_carehome, 
                                             generation_time = generation_time,
                                             estimate_rt = FALSE, fixed = FALSE,
                                             delays = list(incubation_period, reporting_delay_ch),
                                             horizon = 7, samples = 4000, warmup = 500, 
                                             cores = 4, chains = 4, verbose = TRUE, 
                                             adapt_delta = 0.95)

p2 <- estimates_ch$summarised[variable == "infections" & type == "estimate",] %>%
   ggplot2::ggplot(ggplot2::aes(x = date, y= median)) +
   ggplot2::geom_line(col = "dodgerblue") +
   ggplot2::geom_ribbon(data = reported_cases_carehome, ggplot2::aes(x = date, ymax = confirm, ymin = 0), 
                        inherit.aes = FALSE, lty = 2, fill = "red4", alpha = 0.8) +
   ggplot2::geom_ribbon(ggplot2::aes(ymin = bottom, ymax = top), alpha = 0.25, fill = "dodgerblue") +
   ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.25, fill = "dodgerblue") +
   ggplot2::geom_vline(xintercept = as.Date("2020-03-23"), lty = 2) +
   cowplot::theme_cowplot() +
   ggplot2::labs(y = "Daily incidence", x = "Date")

### Joint plot

cm <- estimates$summarised[variable == "infections" & type == "estimate"][, location := "community"]
ch <- estimates_ch$summarised[variable == "infections" & type == "estimate"][ , location := "carehomes"]



df1 <- merge(cm, reported_cases_community[,location := "community"], all = TRUE, by = c("date","location"))
df2 <- merge(ch, reported_cases_carehome[,location := "carehomes"], all = TRUE, by = c("date", "location"))

rbind(df1, df2) %>%
   ggplot2::ggplot(ggplot2::aes(x = date, y= median)) +
   ggplot2::geom_line(col = "dodgerblue") +
   ggplot2::geom_ribbon(ggplot2::aes(x = date, ymax = confirm, ymin = 0), 
                        inherit.aes = FALSE, lty = 2, fill = "red4", alpha = 0.8) +
   ggplot2::geom_ribbon(ggplot2::aes(ymin = bottom, ymax = top), alpha = 0.25, fill = "dodgerblue") +
   ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.25, fill = "dodgerblue") +
   ggplot2::geom_vline(xintercept = as.Date("2020-03-23"), lty = 2) +
   cowplot::theme_cowplot() +
   ggplot2::labs(y = "Daily incidence", x = "Date") +
   ggplot2::facet_wrap( ~ location, ncol = 1)


hi <- merge(cm[, .(date, median)], ch[, .(date, median)], by = "date")[, median := median.x + median.y]

daty <- data.table(date = seq.Date(from = as.Date("2020-05-14"), to = as.Date("2020-07-30"), by = "7 days"),
                   infections = c(4100, 3100, 2500, 2100, 1900, 1800, 1800, 1900, 2000, 2400, 2900, 3700),
                   top = c(6400, 4300, 3300, 2800, 2500, 2400, 2300, 2500, 2700, 3300, 4300, 6400),
                   bottom = c(2500, 2100, 1800, 1500, 1400, 1300, 1300, 1400, 1500, 1600, 1900, 2100))

IFR <- 0.015

cm[, .(date, infections = median / IFR, top = top / IFR, bottom = bottom / IFR)] %>%
   ggplot2::ggplot(ggplot2::aes(x = date, y = infections, ymin = bottom, ymax = top)) +
   ggplot2::geom_ribbon(alpha = 0.5) +
   ggplot2::geom_line() +
   ggplot2::geom_point(data = daty) +
   ggplot2::geom_errorbar(data = daty) +
   ggplot2::labs(y = "Daily new infections", x = "Date") +
   cowplot::theme_cowplot()

pd <- cm[, .(date, infections = median / IFR, top = top / IFR, bottom = bottom / IFR)
][, .(date, infections = cumsum(infections), top = cumsum(top), bottom = cumsum(bottom))]

pd %>%
   ggplot2::ggplot(ggplot2::aes(x = date, y = infections / 56000000,
                                ymin = bottom / 56000000, ymax = top / 56000000)) +
   ggplot2::geom_ribbon(alpha = 0.5) +
   ggplot2::geom_line() +
   # ggplot2::geom_point(data = daty) +
   # ggplot2::geom_errorbar(data = daty) +
   ggplot2::labs(y = "Attack Rate (%)", x = "Date") +
   # ggplot2::ggtitle("IFR = 1.5%") +
   ggplot2::scale_y_continuous(breaks = seq(0, 0.05, 0.01), labels = seq(0:5)) +
   cowplot::theme_cowplot()

