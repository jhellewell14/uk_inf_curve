# Load in libraries
library(data.table)
library(magrittr)
library(ggplot2)
library(scales)
library(patchwork)
library(rstanarm)
options(mc.cores = parallel::detectCores())

# READ ONS LINELIST from covid19_automation
path_to_factory <- "~/repos/covid19_automation"
file_path <- file.path(path_to_factory, "data", "rds", "deaths_eng_latest.rds")
key <- cyphr::data_key(file.path(path_to_factory, "data"))
x <- cyphr::decrypt(readRDS(file_path), key)

ons_linelist <- data.table::as.data.table(x)

# Binary variable for care home or hospital death
ons_linelist[, care_home_death := fifelse(residence_type == "care_nursing_home" |
                                            place_of_death == "care_home",
                                          "Care home",
                                          "Other")]

# Vectors for age groups
agebreaks <- c(0, seq(20, 85, 5), 90, 100)
agelabs <- c("0-19", paste0(seq(20, 85, 5),"-", seq(24, 89, 5)), "90-100")
# Young groups don't have many deaths so they are grouped together for onset to death delays
young_groups <- agelabs[1:4]
old_ind <- 5

# Assign everyone with an age into their age group
ons_linelist[, age_grp := cut(age, breaks = agebreaks, labels = agelabs, right = FALSE)
]

# Assigns missing age groups randomly using probabilities equal to the proportions of people in age groups already
age_dist <- ons_linelist[, .N, age_grp][order(age_grp)]
probs <- age_dist[!is.na(age_grp), N]/sum(age_dist[!is.na(age_grp), N])
ons_linelist[is.na(age_grp), age_grp := sample(x = agelabs, size = age_dist[is.na(age_grp), N], prob = probs, replace = TRUE)]


# Create IFR table stratified by age group
# From here: https://www.medrxiv.org/content/10.1101/2020.07.23.20160895v4.full.pdf
# Takes median age of age group and uses it in meta-regression fit
ifr_tab <- ons_linelist[!is.na(age)]
ifr_tab <- ifr_tab[, .(age = median(age)), age_grp]

ifr_meta <- data.table::fread("data/meta-regression-results.csv")
ifr_meta <- ifr_meta[, .(age, ifr = ifr / 100, ifr_upper = upper / 100, ifr_lower  = lower / 100)]

IFR <- merge(ifr_meta, ifr_tab, by = "age")[order(age_grp)]

## READ CO-CIN LINELIST
# Read in data
# path_to_cocin <- "~/Downloads/CCPUKSARI_DATA_2020-09-14_1105.csv"
path_to_cocin <- "~/Downloads/CCPUKSARI_DATA_2020-08-04_0947.csv"
data <- data.table::fread(path_to_cocin, na.strings = "")

# Select columns + fix read in issue where entries are "" instead of NA
cocin_linelist <- data[,.(cestdat, dsstdtc, dsterm, subjid, age = age_estimateyears)]

cocin_linelist[, c("onset_date_missing", "outcome_date_missing", "dead") :=
                 list(all(is.na(cestdat)), all(is.na(dsstdtc)), any(dsterm == 4, na.rm = TRUE)), 
               by = "subjid"]

# Select people that died with onset dates and death dates
# Calculate number of days between onset and death
# Randomly sample ~ uniform(delay, delay + 1) for delay fit
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

# Assign missing age groups with same probabilities as earlier
cocin_linelist <- cocin_linelist[is.na(age_grp), 
                                 age_grp := sample(agelabs, 
                                                   size = sum(is.na(cocin_linelist$age_grp)), 
                                                   replace = TRUE, prob = probs)][order(age_grp)]


## Fit bootstrapped onset to death delay for young age groups
young_delay <- EpiNow2::bootstrapped_dist_fit(values = cocin_linelist[age_grp %in% young_groups, delay_sampled], 
                              bootstraps = 10,
                              bootstrap_samples = 100,
                              verbose = TRUE)

## Loop over and fit onset to death delay for older age groups
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

## Put together linelist of deaths in the community stratified by age
deaths_community <- ons_linelist[ons == "reported_by_ons" & care_home_death == "Other", 
                          .(confirm = .N, date = date_death), by = c("age_grp", "date_death")
][,.(age_grp, date, confirm)]

## Fill out and assign zero to dates where there were no deaths in an age group
min_date <- min(deaths_community$date) - 10
max_date <- max(deaths_community$date)
deaths_community <- deaths_community[deaths_community[, .(date = seq.Date(from = min_date, to = max_date, by = "day")),
                                 by = .(age_grp)],
                    on = .(age_grp, date),
                    roll = 0][is.na(confirm), confirm := 0][order(age_grp, date)]


## Define generation time and incubation period
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

## Fit EpiNow2 to each age group
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
                                            horizon = 0, samples = 4000, warmup = 500, 
                                            cores = 4, chains = 4, verbose = TRUE, 
                                            adapt_delta = 0.95)
  
  res[[i]] <- estimates$summarised[, age_grp := agelabs[i]][, location := "community"]
  samps[[i]] <- estimates$samples[, age_grp := agelabs[i]][, location := "community"]
}

# Combine results
fr <- data.table::rbindlist(res)

# Sense check plots
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



## Plot all three curves

joint_out <- rbindlist(list(fr[type == "estimate" & variable == "infections"
                               ][, .(median = sum(median), top = sum(top), bottom = sum(bottom)), by = c("date", "location")], 
               fr2[type == "estimate" & variable == "infections", .(date, median, top, bottom, location)]), use.names = TRUE)

## Plot

joint_out %>%
  ggplot(aes(x = date, y = median, ymin = bottom, ymax = top)) +
  geom_line(aes(col = location)) +
  geom_ribbon(aes(fill = location), alpha = 0.5) +
  cowplot::theme_cowplot() +
  labs(x = "Date", y = "Daily infections (that lead to deaths)") +
  geom_vline(xintercept = as.Date("2020-03-23"), lty = 2)

## Join IFR by age group onto EpiNow2 output 

temp_ch <- ons_linelist[ons == "reported_by_ons" & care_home_death == "Care home"]
setkey(temp_ch, age_grp)
temp_ch <- temp_ch[levels(age_grp), .N, by = .EACHI]

IFR$community <- ons_linelist[ons == "reported_by_ons" & care_home_death == "Other"][, .N, by = age_grp][, prop := N / sum(N)][order(age_grp)][, prop]
IFR$carehome <- temp_ch$N / sum(temp_ch$N)

IFR <- melt(IFR, id.vars = c("age_grp", "ifr", "ifr_lower", "ifr_upper"), measure.vars = c("community", "carehome"))[, .(age_grp, location = variable, ifr, ifr_upper, ifr_lower)]
setkey(IFR, age_grp, location)

final_out <- fr[type == "estimate" & variable == "infections"]

# Pad out early zeros in young age groups that EpiNow2 has stripped out
final_out <- final_out[final_out[, .(date = seq.Date(from = min(final_out$date), to = max(final_out$date), by = "day")),
                                 by = .(age_grp)],
                       on = .(age_grp, date),
                       roll = 0][is.na(median), c("median", "top", "bottom", "location") := list(0, 0, 0, "community")]


## Function to smooth youngest age group cases 
smooth_agegrp <- function(x, win = 21){
  y <- c()
  for(i in 1:length(x)) {
    y[i] <- mean(x[(ifelse(i - win <= 0, 1, i - win)):ifelse(i + win > length(x), length(x), i + win)])
  }
  return(y)
}


## Smooth young age groups
final_out[age_grp %in% young_groups, median := smooth_agegrp(median), by = age_grp]
final_out[age_grp %in% young_groups, top := smooth_agegrp(top), by = age_grp]
final_out[age_grp %in% young_groups, bottom := smooth_agegrp(bottom), by = age_grp]


final_out[age_grp %in% young_groups] %>%
  ggplot(aes(x = date, y = median)) +
  geom_line() +
  facet_wrap(~ age_grp, scales = "free_y") +
  geom_ribbon(aes(ymin = bottom, ymax = top), alpha = 0.5)

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
population <- fread(here::here("data/england_population.csv"), header = TRUE)

population <- population[,.(age_grp, population = pop2020, age_upper, age_lower)]

population$age_grp <- as.factor(population$age_grp)

## React 1 and React 2 have slightly different age groups so we need to create 2 
## copies of the population table with merged populations

# React 1
# We can't do 18-24 due to age groups in population data so start at 25-34
pop_react1 <- population[age_lower >= 25]

# Re-factor age groups
pop_react1$age_grp <- rockchalk::combineLevels(pop_react1$age_grp, levs = c("25-29", "30-34"), newLabel = "25-34")
pop_react1$age_grp <- rockchalk::combineLevels(pop_react1$age_grp, levs = c("35-39", "40-44"), newLabel = "35-44")
pop_react1$age_grp <- rockchalk::combineLevels(pop_react1$age_grp, levs = c("45-49", "50-54"), newLabel = "45-54")
pop_react1$age_grp <- rockchalk::combineLevels(pop_react1$age_grp, levs = c("55-59", "60-64"), newLabel = "55-64")
pop_react1$age_grp <- rockchalk::combineLevels(pop_react1$age_grp, levs = c("65-69", "70-74", "75-79", "80-84", "85-89", "90+"), newLabel = "65-100")

# Sum up by new age groups
pop_react1 <- pop_react1[, .(age = sum(population)), by = age_grp]

# React 2
# Again, we can't do 18-24 due to age groups in population data so start at 25-34
pop_react2 <- population[age_lower >= 25]

# Re-factor age groups
pop_react2$age_grp <- rockchalk::combineLevels(pop_react2$age_grp, levs = c("25-29", "30-34"), newLabel = "25-34")
pop_react2$age_grp <- rockchalk::combineLevels(pop_react2$age_grp, levs = c("35-39", "40-44"), newLabel = "35-44")
pop_react2$age_grp <- rockchalk::combineLevels(pop_react2$age_grp, levs = c("45-49", "50-54"), newLabel = "45-54")
pop_react2$age_grp <- rockchalk::combineLevels(pop_react2$age_grp, levs = c("55-59", "60-64"), newLabel = "55-64")
pop_react2$age_grp <- rockchalk::combineLevels(pop_react2$age_grp, levs = c("65-69", "70-74"), newLabel = "65-74")
pop_react2$age_grp <- rockchalk::combineLevels(pop_react2$age_grp, levs = c("75-79", "80-84", "85-89", "90+"), newLabel = "75-100")

# Sum up by new age groups
pop_react2 <- pop_react2[, .(population = sum(population)), by = age_grp]

## RE-FORMULATE EPINOW2 OUTPUT TO REACT 1 & 2 AGE GROUPS

# React 1
final_out_react1 <- final_out[, .(age_grp, date, median = median / ifr, top = top / ifr_lower, bottom = bottom / ifr_upper)
                              ][!(age_grp %in% c("0-19", "20-24"))]
# Re-factor age groups
final_out_react1$age_grp <- as.factor(final_out_react1$age_grp)
final_out_react1$age_grp <- rockchalk::combineLevels(final_out_react1$age_grp, levs = c("25-29", "30-34"), newLabel = "25-34")
final_out_react1$age_grp <- rockchalk::combineLevels(final_out_react1$age_grp, levs = c("35-39", "40-44"), newLabel = "35-44")
final_out_react1$age_grp <- rockchalk::combineLevels(final_out_react1$age_grp, levs = c("45-49", "50-54"), newLabel = "45-54")
final_out_react1$age_grp <- rockchalk::combineLevels(final_out_react1$age_grp, levs = c("55-59", "60-64"), newLabel = "55-64")
final_out_react1$age_grp <- rockchalk::combineLevels(final_out_react1$age_grp, levs = c("65-69", "70-74", "75-79", "80-84", "85-89", "90-100"), newLabel = "65-100")

# Sum back up
final_out_react1 <- final_out_react1[, .(age_grp, date, bottom, top, median)][, .(bottom = sum(bottom), median = sum(median), top = sum(top)), by = c("age_grp", "date")]

# React 2
final_out_react2 <- final_out[, .(age_grp, date, median = median / ifr, top = top / ifr_lower, bottom = bottom / ifr_upper)
                              ][!(age_grp %in% c("0-19", "20-24"))]
# Re-factor age groups
final_out_react2$age_grp <- as.factor(final_out_react2$age_grp)
final_out_react2$age_grp <- rockchalk::combineLevels(final_out_react2$age_grp, levs = c("25-29", "30-34"), newLabel = "25-34")
final_out_react2$age_grp <- rockchalk::combineLevels(final_out_react2$age_grp, levs = c("35-39", "40-44"), newLabel = "35-44")
final_out_react2$age_grp <- rockchalk::combineLevels(final_out_react2$age_grp, levs = c("45-49", "50-54"), newLabel = "45-54")
final_out_react2$age_grp <- rockchalk::combineLevels(final_out_react2$age_grp, levs = c("55-59", "60-64"), newLabel = "55-64")
final_out_react2$age_grp <- rockchalk::combineLevels(final_out_react2$age_grp, levs = c("65-69", "70-74"), newLabel = "65-74")
final_out_react2$age_grp <- rockchalk::combineLevels(final_out_react2$age_grp, levs = c("75-79", "80-84", "85-89", "90-100"), newLabel = "75-100")

# Sum back up
final_out_react2 <- final_out_react2[, .(age_grp, date, bottom, top, median)][, .(bottom = sum(bottom), median = sum(median), top = sum(top)), by = c("age_grp", "date")]

### VERY SIMPLE COMPARTMENTAL DECAY APPROACH
decay_inf <- function(x, decay_rate, test_sens, test_spec) {
  out <- c()
  out[1] <- x[1]
  for(t in 2:length(x)) {
    out[t] <- out[t - 1] + x[t] - (decay_rate * out[t - 1])
  }
  return(out * test_sens)
}

# Plot REACT 1 results
# PCR test parameters
av_test_neg <- 5 # this is value PHE use via Nick
pcr_sensitivity <- 1
pcr_specificity <- 1

final_out_react1 <- merge(final_out_react1, pop_react1, by = "age_grp")[order(age_grp, date)]

# Create decayed prevalence variables
final_out_react1[, dec_prev := decay_inf(median, decay_rate = 1 / av_test_neg, test_sens = pcr_sensitivity, test_spec = pcr_specificity), by = age_grp]
final_out_react1[, dec_bot := decay_inf(bottom, decay_rate = 1 / av_test_neg, test_sens = pcr_sensitivity, test_spec = pcr_specificity), by = age_grp]
final_out_react1[, dec_top := decay_inf(top, decay_rate = 1 / av_test_neg, test_sens = pcr_sensitivity, test_spec = pcr_specificity), by = age_grp]

final_out_react1 %>%
  ggplot(aes(x = date, y = dec_prev / age, ymin = dec_bot / age, ymax = dec_top / age)) +
  geom_line() +
  geom_ribbon(alpha = 0.4) +
  # scale_y_continuous(labels = comma) +
  # geom_vline(xintercept = as.Date("2020-03-23")) +
  cowplot::theme_minimal_grid() +
  ggplot2::labs(x = "Date", y = "Prevalence by swab + PCR (%)", title = paste0("Average time until PCR-negative: ", av_test_neg, " days after infection")) +
  facet_wrap(~ age_grp) +
  geom_errorbarh(data = subset(sero, study == "React 1" & age_grp != "18-24"),
                 inherit.aes = FALSE, aes(xmin = as.Date(start_date), xmax = as.Date(end_date),
                                          y = seroprev / 100, group = age_grp), col = "red4") +
  geom_errorbar(data = subset(sero, study == "React 1" & age_grp != "18-24"),
                inherit.aes = FALSE, aes(x = start_date + (end_date - start_date) / 2, ymin = lower / 100, ymax = upper / 100),
                col = "red4", width = 0) +
  scale_y_continuous(breaks = seq(0, 0.04, 0.01), labels = seq(0, 4, 1)) +
  coord_cartesian(xlim = c(as.Date("2020-05-01"), as.Date("2020-09-07")))


## ALTERNATIVE REACT 1 PLOT
plot_dt <- merge(final_out_react1[age_grp != "0-34"], sero[study == "React 1" & age_lower >= 35 & round <= 3], by = "age_grp", allow.cartesian = TRUE)
plot_dt[date >= start_date & date <= end_date
][, .(median = median(dec_prev / age), 
      top = median(dec_top / age),
      bottom = median(dec_bot / age),
      prev = unique(seroprev) / 100, 
      lower = unique(lower) / 100, 
      upper = unique(upper) / 100), list(age_grp, round)] %>%
  ggplot() +
  geom_errorbar(aes(x = age_grp, ymin = bottom, ymax = top)) + 
  geom_errorbar(aes(x = age_grp, ymin = lower, ymax = upper), col = "red4") +
  facet_wrap(~ round) +
  cowplot::theme_minimal_grid() +
  scale_y_continuous(breaks = seq(0, 0.005, 0.001), labels = seq(0, 0.5, 0.1)) +
  labs(y = "PCR prevalence (%)", x = "Age group")

## PLOT RESULTS FOR REACT 2

# Average time to sero-reversion
av_sero <- 180 / log(2) # calculating mean time to sero-reversion from 180 day half life
sero_sensitivity <- 1
sero_specificity <- 1

final_out_react2 <- merge(final_out_react2, pop_react2, by = "age_grp")[order(age_grp, date)]

final_out_react2[, dec_inf := decay_inf(median, 1 / av_sero, 1, 1), age_grp]
final_out_react2[, dec_bot := decay_inf(bottom, 1 / av_sero, 1, 1), age_grp]
final_out_react2[, dec_top := decay_inf(top, 1 / av_sero, 1, 1), age_grp] 

final_out_react2 %>%
  ggplot(aes(x = date, y = dec_inf / population, ymin = dec_bot / population, ymax = dec_top / population)) + 
  geom_line() + 
  geom_ribbon(alpha = 0.4) + 
  facet_wrap(~ age_grp) +
  geom_point(data = subset(sero, study == "React 2" & age_lower >= 25), aes(x = start_date + (end_date - start_date) / 2, y = seroprev / 100), col = "red4", inherit.aes = FALSE) +
  geom_errorbar(data = subset(sero, study == "React 2" & age_lower >= 25), aes(x = start_date + (end_date - start_date) / 2, ymin = lower / 100, ymax = upper / 100), col = "red4", inherit.aes = FALSE) +
  geom_errorbarh(data = subset(sero, study == "React 2" & age_lower >= 25), aes(xmin = start_date, xmax = end_date, y = seroprev / 100), col = "red4", inherit.aes = FALSE) +
  scale_y_continuous(breaks = seq(0, 0.13, 0.01), labels = seq(0, 13, 1)) +
  labs(y = "Seroprevalence (%)", x = "Date", title  = paste0("Average time to sero-reversion: ", av_sero, " days")) +
  cowplot::theme_minimal_grid()




