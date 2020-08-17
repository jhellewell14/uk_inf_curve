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

ons_linelist[is.na(age_grp), age_grp := "80-89"]


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



