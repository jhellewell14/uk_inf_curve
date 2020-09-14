library(data.table)
library(magrittr)
library(ggplot2)
library(scales)
library(patchwork)

# READ ONS LINELIST
path_to_factory <- "~/repos/covid19_automation"
file_path <- file.path(path_to_factory, "data", "rds", "deaths_eng_latest.rds")
key <- cyphr::data_key(file.path(path_to_factory, "data"))
x <- cyphr::decrypt(readRDS(file_path), key)

ons_linelist <- data.table::as.data.table(x)

ons_linelist[, care_home_death := fifelse(residence_type == "care_nursing_home" |
                                            place_of_death == "care_home",
                                          "Care home",
                                          "Other")]


agebreaks <- c(0, 35, 45, 55, 65, 75, 100)
agelabs <- c("0-34", "35-44", "45-54", "55-64", "65-74","75-100")
young_groups <- "0-34"
old_ind <- 2

midpoints <- c()
for(i in 2:length(agebreaks)){
  midpoints[i - 1] <- (agebreaks[i - 1] + agebreaks[i]) / 2
}

source(here::here("R/fit_IFR_agegroup.R"))

preds <- predict(fit, newdata = data.frame(agemid = midpoints), se.fit = TRUE)

IFR <- data.table(age_grp = agelabs,
           ifr = boot::inv.logit(preds$fit),
           ifr_lower = boot::inv.logit(preds$fit - 1.96*preds$se.fit),
           ifr_upper = boot::inv.logit(preds$fit + 1.96*preds$se.fit))


ons_linelist[, age_grp := cut(age, breaks = agebreaks, labels = agelabs, right = FALSE)
]

# Assigns missing age groups randomly using age distribution found in the data set
# Maybe a better way of doing this
age_dist <- ons_linelist[, .N, age_grp][order(age_grp)]
probs <- age_dist[!is.na(age_grp), N]/sum(age_dist[!is.na(age_grp), N])
ons_linelist[is.na(age_grp), age_grp := sample(x = agelabs, size = age_dist[is.na(age_grp), N], prob = probs, replace = TRUE)]

## READ CO-CIN LINELIST
# Read in data
data <- data.table::fread("~/Downloads/CCPUKSARI_DATA_2020-08-04_0947.csv", na.strings = "")

# Select columns + fix read in issue where entries are "" instead of NA
cocin_linelist <- data[,.(cestdat, dsstdtc, dsterm, subjid, age = age_estimateyears)]

cocin_linelist[, c("onset_date_missing", "outcome_date_missing", "dead") :=
                 list(all(is.na(cestdat)), all(is.na(dsstdtc)), any(dsterm == 4, na.rm = TRUE)), 
               by = "subjid"]

cocin_linelist <- cocin_linelist[!onset_date_missing & !outcome_date_missing & dead
                                 ][, .(onset_date = unique(cestdat[!is.na(cestdat)]),
                                       dead = unique(dead),
                                       age = unique(age[!is.na(age)]),
                                       outcome_date = unique(dsstdtc[!is.na(dsstdtc)])), by = "subjid"
                                   ][, delay := as.integer(as.Date(outcome_date) - as.Date(onset_date))
                                     ][delay >= 0 & delay <= 60 & as.Date(onset_date) > "2020-01-01"
                                       ][, delay_sampled := runif(.N, delay, delay + 1)]

# This again assigns NA age groups by distribution in data again
cocin_linelist <- cocin_linelist[, age_grp := cut(age, breaks = agebreaks, labels = agelabs, right = FALSE)]

cocin_linelist <- cocin_linelist[is.na(age_grp), age_grp := sample(agelabs, size = sum(is.na(cocin_linelist$age_grp)), replace = TRUE, prob = probs)][order(age_grp)]

## SAMPLE REPORTING DELAYS BY AGE

young_delay <- EpiNow2::bootstrapped_dist_fit(values = cocin_linelist[age_grp %in% young_groups, delay_sampled], 
                              bootstraps = 10,
                              bootstrap_samples = 100,
                              verbose = TRUE)

delays <- list()
for(i in 1:length(agelabs)) {
  print(i)
  if(i < old_ind){
    delays[[i]] <- young_delay
    delays[[i]]$age_grp <- young_groups[i]
  }else{
    delays[[i]] <- EpiNow2::bootstrapped_dist_fit(values = cocin_linelist[age_grp == agelabs[i], delay_sampled],
                                                  bootstraps = 10,
                                                  bootstrap_samples = 100,
                                                  verbose = TRUE)
    delays[[i]]$age_grp <- agelabs[i]
  }
}
## COMMUNITY DEATHS BY AGE
deaths_community <- ons_linelist[ons == "reported_by_ons" & care_home_death == "Other", 
                          .(confirm = .N, date = date_death), by = c("age_grp", "date_death")
][,.(age_grp, date, confirm)]

deaths_community <- deaths_community[deaths_community[, .(date = seq.Date(from = min(date), to = max(date), by = "day")),
                                 by = .(age_grp)],
                    on = .(age_grp, date),
                    roll = 0][is.na(confirm), confirm := 0][order(age_grp, date)]


## Fit EpiNow2 by age
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

res <- list()
samps <- list()
for(i in 1:length(agelabs)){
  print(i)
  
  reporting_delay <- as.list(rbindlist(delays)[age_grp == agelabs[i]])
  reporting_delay$max <- 60
  
  estimates <- EpiNow2::estimate_infections(reported_cases = deaths_community[age_grp == agelabs[i]], 
                                            generation_time = generation_time,
                                            estimate_rt = FALSE, fixed = FALSE,
                                            delays = list(incubation_period, reporting_delay),
                                            horizon = 7, samples = 4000, warmup = 500, 
                                            cores = 4, chains = 4, verbose = TRUE, 
                                            adapt_delta = 0.95)
  
  res[[i]] <- estimates$summarised[, age_grp := agelabs[i]][, location := "community"]
  samps[[i]] <- estimates$samples[, age_grp := agelabs[i]][, location := "community"]
}

fr <- data.table::rbindlist(res)


fr[type == "estimate" & variable == "infections"] %>%
  ggplot2::ggplot(ggplot2::aes(x = date, y = median, col = age_grp)) +
  ggplot2::geom_line() +
  ggplot2::scale_color_discrete(name = "Age group") +
  cowplot::theme_cowplot() +
  ggplot2::labs(x = "Date", y = "Infections that lead to deaths")

fr[type == "estimate" & variable == "infections"][, .(median = sum(median), top = sum(top), bottom = sum(bottom)), by = "date"] %>%
  ggplot2::ggplot(ggplot2::aes(x = date, y = median, ymin = bottom, ymax = top)) +
  ggplot2::geom_line() +
  ggplot2::geom_ribbon(alpha = 0.5) +
  cowplot::theme_cowplot() +
  ggplot2::labs(x = "Date", y = "Infections that lead to deaths")

## CARE HOME DEATHS BY AGE

deaths_carehome <- ons_linelist[ons == "reported_by_ons" & care_home_death == "Care home", 
                                .(confirm = .N, date = date_death), by = c("age_grp", "date_death")
][,.(age_grp, date, confirm)]


deaths_carehome <- deaths_carehome[deaths_carehome[, .(date = seq.Date(from = min(date), to = max(date), by = "day")),
                                                      by = .(age_grp)],
                                     on = .(age_grp, date),
                                     roll = 0][is.na(confirm), confirm := 0][order(age_grp, date)]

carehome_delays <- ons_linelist[ons == "reported_by_ons" & care_home_death == "Care home" & !is.na(date_death) & !is.na(date_onset)
             ][, delay := as.integer(as.Date(date_death) - as.Date(date_onset))
               ][delay >= 0 & delay <= 60 & as.Date(date_onset) > "2020-01-01"
               ][, delay_sampled := runif(.N, delay, delay + 1)
                 ][, .(delay_sampled, age_grp)][order(age_grp)]

carehome_delays <- EpiNow2::bootstrapped_dist_fit(values = carehome_delays[, delay_sampled],
                                                  bootstraps = 10,
                                                  bootstrap_samples = 100,
                                                  verbose = TRUE)

carehome_fit <- estimates <- EpiNow2::estimate_infections(reported_cases = deaths_carehome[, .(confirm = sum(confirm)), by = "date"], 
                                                            generation_time = generation_time,
                                                            estimate_rt = FALSE, fixed = FALSE,
                                                            delays = list(incubation_period, carehome_delays),
                                                            horizon = 7, samples = 4000, warmup = 500, 
                                                            cores = 4, chains = 4, verbose = TRUE, 
                                                            adapt_delta = 0.95)

fr2 <- carehome_fit$summarised[, location := "Care home"]

## Hospital acquired infection

hosp_acq_prop <- 0.11

hosp <- fr[type == "estimate" & variable == "infections" & location == "community"]

hosp[, location := "hospital"]

cols <- c("bottom", "top", "lower", "upper", "median", "mean")
hosp[, (cols) := lapply(.SD, "*", hosp_acq_prop), .SDcols = cols]


fr[type == "estimate" & variable == "infections" & location == "community", (cols) := lapply(.SD, "*", 1 - hosp_acq_prop), .SDcols = cols]

fr <- rbind(fr, hosp)

## Plot all three curves

joint_out <- rbindlist(list(fr[type == "estimate" & variable == "infections"
                               ][, .(median = sum(median), top = sum(top), bottom = sum(bottom)), by = c("date", "location")], 
               fr2[type == "estimate" & variable == "infections", .(date, median, top, bottom, location)]), use.names = TRUE)

joint_out %>%
  ggplot(aes(x = date, y = median, ymin = bottom, ymax = top)) +
  geom_line(aes(col = location)) +
  geom_ribbon(aes(fill = location), alpha = 0.5) +
  cowplot::theme_cowplot() +
  labs(x = "Date", y = "Daily infections (that lead to deaths)") +
  geom_vline(xintercept = as.Date("2020-03-23"), lty = 2)

## IFR 

temp_ch <- ons_linelist[ons == "reported_by_ons" & care_home_death == "Care home"]
setkey(temp_ch, age_grp)
temp_ch <- temp_ch[levels(age_grp), .N, by = .EACHI]

IFR$community <- ons_linelist[ons == "reported_by_ons" & care_home_death == "Other"][, .N, by = age_grp][, prop := N / sum(N)][order(age_grp)][, prop]
IFR$carehome <- temp_ch$N / sum(temp_ch$N)

IFR <- melt(IFR, id.vars = c("age_grp", "ifr", "ifr_lower", "ifr_upper"), measure.vars = c("community", "carehome"))[, .(age_grp, location = variable, ifr, ifr_upper, ifr_lower)]
setkey(IFR, age_grp, location)

final_out <- fr[type == "estimate" & variable == "infections"]


## Smooth youngest age group cases
x <- final_out[age_grp == "0-39", median]

y <- c()
win <- 7
for(i in 1:length(x)) {
  y[i] <- mean(x[(ifelse(i - win <= 0, 1, i - win)):ifelse(i + win > length(x), length(x), i + win)])
}

final_out[age_grp == "0-39", median := y]

setkey(final_out, date, age_grp, location)

## Final plots
final_out <- merge(final_out, IFR, by = c("age_grp", "location"))

p1 <- final_out[, .(date, median = median / ifr, top = top / ifr_lower, bottom = bottom / ifr_upper)
          ][, .(median = sum(median), top = sum(top), bottom = sum(bottom)), by = "date"] %>%
  ggplot(aes(x = date, y = median, ymin = bottom, ymax = top)) +
  geom_line() +
  geom_ribbon(alpha = 0.4) +
  scale_y_continuous(labels = comma) +
  geom_vline(xintercept = as.Date("2020-03-23")) +
  cowplot::theme_cowplot() +
  ggplot2::labs(x = "Date", y = "Daily infections")

p2 <- final_out[, .(date, median = median / ifr, top = top / ifr_lower, bottom = bottom / ifr_upper)
][, .(median = sum(median), top = sum(top), bottom = sum(bottom)), by = "date"
  ][order(date)][, .(median = cumsum(median), top = cumsum(top), bottom = cumsum(bottom), date)] %>%
  ggplot(aes(x = date, y = median, ymin = bottom, ymax = top)) +
  geom_line() +
  geom_ribbon(alpha = 0.4) + 
  scale_y_continuous(labels = comma) + 
  geom_vline(xintercept = as.Date("2020-03-23")) +
  cowplot::theme_cowplot() +
  ggplot2::labs(x = "Date", y = "Cumulative daily infections")

p1 + p2


## READ IN AND FORMAT REACT STUDY DATA

sero <- fread(here::here("data/react_results.csv"))[study %in% c("React 1", "React 2")]

sero <- sero[, age_grp := paste0(age_lower, "-", age_upper)
             ][age_upper >= 18]

sero[, start_date := as.Date(start_date, format = "%d-%m-%Y")
     ][, end_date := as.Date(end_date, format = "%d-%m-%Y")]



## READ IN POPULATION DATA
population <- fread(here::here("data/age_collated.csv"), header = TRUE)

population$age_grp <- as.factor(population$age_grp)


final_out <- merge(final_out, population, by = "age_grp")


### VERY SIMPLE COMPARTMENTAL DECAY APPROACH
decay_inf <- function(x, decay_rate, test_sens, test_spec) {
  out <- c()
  out[1] <- x[1]
  for(t in 2:length(x)) {
    out[t] <- out[t - 1] + x[t] - (decay_rate * out[t - 1])
  }
  return(out * test_sens)
}


## PLOT RESULTS FOR REACT 1 (SWAB) STUDY
## This requires merging age groups
pop_react1 <- population
pop_react1$age_grp <- rockchalk::combineLevels(pop_react1$age_grp, levs = c("65-74", "75-100"), newLabel = "65-100")
pop_react1 <- pop_react1[, .(age = sum(Age_2020)), by = age_grp]

final_out_react1 <- final_out[, .(age_grp, date, median = median / ifr, top = top / ifr_lower, bottom = bottom / ifr_upper)]

final_out_react1$age_grp <- as.factor(final_out_react1$age_grp)
final_out_react1$age_grp <- rockchalk::combineLevels(final_out_react1$age_grp, levs = c("65-74", "75-100"), newLabel = "65-100")

final_out_react1 <- final_out_react1[, .(age_grp, date, bottom, top, median)][, .(bottom = sum(bottom), median = sum(median), top = sum(top)), by = c("age_grp", "date")]

final_out_react1 <- merge(final_out_react1, pop_react1, by = "age_grp")

# Average time after infection to test negative
av_test_neg <- 10
pcr_sensitivity <- 0.83
pcr_specificity <- 0.933

final_out_react1[, dec_prev := decay_inf(median, decay_rate = 1 / av_test_neg, test_sens = pcr_sensitivity, test_spec = pcr_specificity), by = age_grp]
final_out_react1[, dec_bot := decay_inf(bottom, decay_rate = 1 / av_test_neg, test_sens = pcr_sensitivity, test_spec = pcr_specificity), by = age_grp]
final_out_react1[, dec_top := decay_inf(top, decay_rate = 1 / av_test_neg, test_sens = pcr_sensitivity, test_spec = pcr_specificity), by = age_grp]

final_out_react1[age_grp != "0-34"] %>%
  ggplot(aes(x = date, y = dec_prev / age, ymin = dec_bot / age, ymax = dec_top / age)) +
  geom_line() +
  geom_ribbon(alpha = 0.4) +
  scale_y_continuous(labels = comma) +
  # geom_vline(xintercept = as.Date("2020-03-23")) +
  cowplot::theme_minimal_grid() +
  ggplot2::labs(x = "Date", y = "Prevalence by swab + PCR (%)", title = paste0("Average time until PCR-negative: ", av_test_neg, " days after infection")) +
  facet_wrap(~ age_grp) +
  geom_errorbarh(data = subset(sero, study == "React 1" & age_grp != "18-24" & age_grp != "25-34"),
                 inherit.aes = FALSE, aes(xmin = as.Date(start_date), xmax = as.Date(end_date), width = 0,
                                          y = seroprev / 100, group = age_grp), col = "red4") +
  geom_errorbar(data = subset(sero, study == "React 1" & age_grp != "18-24" & age_grp != "25-34"),
                inherit.aes = FALSE, aes(x = start_date + (end_date - start_date) / 2, ymin = lower / 100, ymax = upper / 100),
                col = "red4", width = 0) +
  scale_y_continuous(breaks = seq(0, 0.04, 0.01), labels = seq(0, 4, 1))

## PLOT RESULTS FOR REACT 2

# Average time to sero-reversion
av_sero <- 30*6
sero_sensitivity <- 0.9
sero_specificity <- 1

final_out[age_grp != "0-34"
          ][, dec_inf := decay_inf(median / ifr, 1 / av_sero, 1, 1), age_grp
            ][, dec_bot := decay_inf(bottom / ifr_upper, 1 / av_sero, 1, 1), age_grp
              ][, dec_top := decay_inf(top / ifr_lower, 1 / av_sero, 1, 1), age_grp] %>%
  ggplot(aes(x = date, y = dec_inf / Age_2020, ymin = dec_bot / Age_2020, ymax = dec_top / Age_2020)) + 
  geom_line() + 
  geom_ribbon(alpha = 0.4) +
  facet_wrap(~ age_grp) +
  geom_point(data = subset(sero, study == "React 2" & age_lower >= 35), aes(x = start_date, y = seroprev / 100), col = "red4", inherit.aes = FALSE) +
  geom_point(data = subset(sero, study == "React 2" & age_lower >= 35), aes(x = end_date, y = seroprev / 100), col = "red4", inherit.aes = FALSE) +
  geom_errorbar(data = subset(sero, study == "React 2" & age_lower >= 35), aes(x = start_date, ymin = lower / 100, ymax = upper / 100), col = "red4", inherit.aes = FALSE) +
  geom_errorbar(data = subset(sero, study == "React 2" & age_lower >= 35), aes(x = end_date, ymin = lower / 100, ymax = upper / 100), col = "red4", inherit.aes = FALSE) +
  scale_y_continuous(breaks = seq(0, 0.13, 0.01), labels = seq(0, 13, 1)) +
  labs(y = "Seroprevalence (%)", x = "Date", title  = paste0("Average time to sero-reversion: ", av_sero, " days")) +
  cowplot::theme_minimal_grid()




