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

agebreaks <- c(0, 10, 20, 30, 40, 50, 60, 70 ,80, 90, 100)
agelabs <- c("0-9", "10-19", 
             "20-29", "30-39", 
             "40-49", "50-59", 
             "60-69", "70-79",
             "80-89","90-99")

ons_linelist[, age_grp := cut(age, breaks = agebreaks, labels = agelabs, right = FALSE)
]

# Assigns deaths to most likely age group
# Maybe a better way of doing this
ons_linelist[is.na(age_grp), age_grp := "80-89"]

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

cocin_linelist <- cocin_linelist[, age_grp := cut(age, breaks = agebreaks, labels = agelabs, right = FALSE)
                                 ][!is.na(age_grp)][order(age_grp)]

## SAMPLE REPORTING DELAYS BY AGE

young_delay <- EpiNow2::bootstrapped_dist_fit(values = cocin_linelist[age_grp %in% c("0-9", "10-19",
                              "20-29", "30-39"), delay_sampled], 
                              bootstraps = 10,
                              bootstrap_samples = 100,
                              verbose = TRUE)

delays <- list(young_delay, young_delay, young_delay, young_delay)
delays[[1]]$age_grp <- "0-9"
delays[[2]]$age_grp <- "10-19"
delays[[3]]$age_grp <- "20-29"
delays[[4]]$age_grp <- "30-39"

for(i in 1:length(agelabs[5:10])) {
  print(i)
  delays[[i + 4]] <- EpiNow2::bootstrapped_dist_fit(values = cocin_linelist[age_grp == agelabs[4 + i], delay_sampled],
                                              bootstraps = 10,
                                              bootstrap_samples = 100,
                                              verbose = TRUE)
  delays[[i + 4]]$age_grp <- agelabs[4 + i]
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
  
  res[[i]] <- estimates$summarised[, age_grp := agelabs[i]]
  samps[[i]] <- estimates$samples[, age_grp := agelabs[i]]
}

fr <- data.table::rbindlist(res)

fr[type == "estimate" & variable == "infections"] %>%
  ggplot2::ggplot(ggplot2::aes(x = date, y = median, col = age_grp)) +
  ggplot2::geom_line() +
  ggplot2::scale_color_discrete(name = "Age group") +
  cowplot::theme_cowplot() +
  ggplot2::labs(x = "Date", y = "Infections that lead to deaths")



