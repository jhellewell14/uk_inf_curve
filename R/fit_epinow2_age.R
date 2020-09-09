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

## GENEVA ESTIMATE AGE GROUPS AND IFR
# agebreaks <-  c(0, 10, 20, 50, 65, 100)
# agelabs <- c("5-9", "10-19", "20-49", "50-64", "65-100")
# na_fill_group <- "65-100"
# young_groups <- c("5-9", "10-19", "20-49")
# old_ind <- 4
# 
# IFR <- data.table(age_grp = agelabs,
#                   ifr = c(0.000016, 0.0000032, 0.000092, 0.0014, 0.056))

## VERITY ESTIMATE AGE GROUPS AND IFR
# agebreaks <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 100)
# agelabs <- c("0-9", "10-19",
#                           "20-29", "30-39",
#                           "40-49", "50-59",
#                           "60-69", "70-79",
#                           "80-100")
# na_fill_group <- "80-100"
# young_groups <- c("0-9", "10-19", "20-29", "30-39")
# old_ind <- 5
# 
# IFR <- data.table(age_grp = agelabs,
#                   ifr = c(0.00161, 0.00695, 0.0309, 0.0844, 0.161, 0.595, 1.93, 4.28, 7.80) / 100)

## SALJE ESTIMATE AGE GROUPS AND IFR
# agebreaks <- c(0, 20, 30, 40, 50, 60, 70, 80, 100)
# agelabs <- c("0-20", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-100")
# na_fill_group <- "80-100"
# young_groups <- c("0-20", "20-29", "30-39")
# old_ind <- 4
# 
# IFR <- data.table(age_grp = agelabs,
#                   ifr = c(0.001, 0.005, 0.02, 0.05, 0.2, 0.7, 1.9, 8.3) / 100)

# 10 year age groups
# agebreaks <- c(0, 10, 20, 30, 40, 50, 60, 70 ,80, 90, 100)
# agelabs <- c("0-9", "10-19", 
#              "20-29", "30-39", 
#              "40-49", "50-59", 
#              "60-69", "70-79",
#              "80-89","90-99")
# na_fill_group <- "80-89"
# young_groups <- c("0-9", "10-19",
#                   "20-29", "30-39")

# IFR <- data.table(age_grp = agelabs, 
#            # ifr = c(1.6e-05, 7e-05, 0.00031, 0.00084, 0.0016, 0.006, 0.019, 0.043, 0.078, 0.078),
#            ifr = (exp(-8.1290  + (c(5,15,25,35,45,55,65,75,85,95) * 0.1191))) / 100)

## MERGED AGE GROUPS TO PREVENT LOW DEATH PROBLEMS
## IFR FROM LOGIT-LINEAR MODEL
agebreaks <- c(0, 40, 50, 60, 70, 85, 100)
agelabs <- c("0-39", "40-49", "50-59", "60-69", "70-84", "85-100")
young_groups <- "0-39"
old_ind <- 2

midpoints <- c()
for(i in 2:length(agebreaks)){
  midpoints[i - 1] <- (agebreaks[i - 1] + agebreaks[i]) / 2
}

source(here::here("R/fit_IFR_agegroup.R"))
IFR <- data.table(age_grp = agelabs,
           ifr = boot::inv.logit(predict(fit, newdata = data.frame(agemid = midpoints))))


ons_linelist[, age_grp := cut(age, breaks = agebreaks, labels = agelabs, right = FALSE)
]

# Assigns missing age groups randomly using age distribution found in the data set
# Maybe a better way of doing this
age_dist <- ons_linelist[, .N, age_grp]
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

joint_out <- rbindlist(list(fr[type == "estimate" & variable == "infections"
                               ][, .(median = sum(median), top = sum(top), bottom = sum(bottom), location = location[1]), by = date], 
               fr2[type == "estimate" & variable == "infections", .(date, median, top, bottom, location)]))

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

IFR <- melt(IFR, id.vars = c("age_grp", "ifr"), measure.vars = c("community", "carehome"))[, .(age_grp, location = variable, ifr)]
setkey(IFR, age_grp, location)

final_out <- fr[type == "estimate" & variable == "infections"]

x <- final_out[age_grp == "0-39", median]

y <- c()
win <- 7
for(i in 1:length(x)) {
  y[i] <- mean(x[(ifelse(i - win <= 0, 1, i - win)):ifelse(i + win > length(x), length(x), i + win)])
}

final_out[age_grp == "0-39", median := y]

setkey(final_out, date, age_grp, location)

final_out <- merge(final_out, IFR, by = c("age_grp", "location"))

p1 <- final_out[, .(date, median = median / ifr, top = top / ifr, bottom = bottom / ifr)
          ][, .(median = sum(median), top = sum(top), bottom = sum(bottom)), by = "date"] %>%
  ggplot(aes(x = date, y = median, ymin = bottom, ymax = top)) +
  geom_line() +
  geom_ribbon(alpha = 0.4) +
  scale_y_continuous(labels = comma) +
  geom_vline(xintercept = as.Date("2020-03-23")) +
  cowplot::theme_cowplot() +
  ggplot2::labs(x = "Date", y = "Daily infections")

p2 <- final_out[, .(date, median = median / ifr, top = top / ifr, bottom = bottom / ifr)
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

### VERY SIMPLE VALIDATION APPROACH
infections <- final_out[age_grp == "50-59"][, median := median / ifr]
pos_length <- 1 / 91.25 # rate of people becoming sero-negative again
rate_to_neg <- 1 / 91.25 # rate of people beginnning to test PCR-negative
sp <- c()
sp[1] <- 0
sn <- c()
sn[1] <- 0
tp <- c()
tp[1] <- 0
for(t in 2:length(infections$median)) {
  sp[t] <- sp[t - 1] + infections$median[t] - (pos_length * sp[t - 1])
  sn[t] <- sn[t - 1] + (pos_length * sp[t - 1])
  tp[t] <- tp[t - 1] + infections$median[t] - (rate_to_neg * tp[t - 1])
}
plot(infections$date, cumsum(infections$median) / 7578112, type = "l", 
     ylab = "Sero (red) vs PCR (blue) vs True (black) Prevalence",
     xlab = "Date")
lines(infections$date, sp / 7578112, col = "red")
lines(infections$date, tp / 7578112, col = "blue")
