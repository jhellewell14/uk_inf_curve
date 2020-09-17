# Age distribution of hospital-acquired infections

cocin <- data[, .(onset_date = as.Date(cestdat), admission_date = as.Date(hostdat), age = age_estimateyears, outcome = dsterm)
][!is.na(onset_date)
][, hosp_acq := onset_date > (as.Date(admission_date) + 5)]

hosp_age_prop <- cocin[, age_grp := cut(age, breaks = agebreaks, labels = agelabs, right = FALSE)
][!is.na(age_grp) & !is.na(hosp_acq)][, .(y = .N), by = list(hosp_acq, age_grp)
][, .(prop_of_hosp_acq = y / sum(y), age_grp), by = hosp_acq][hosp_acq == TRUE][,hosp_acq := NULL]

# Percentage of infections that require hospitalistion by age

hosp <- fread("data/hospitalisation_rate_verity.csv", header = TRUE)

hosp[, age_mid := (age_higher + age_lower) / 2][age_lower == 0, hosp_rate := 1E-09]

hosp_fit <- lm(data = hosp, hosp_rate ~ age_mid)

preds <- predict(hosp_fit, newdata = data.frame(age_mid = seq(17, 90, 1)), se.fit = TRUE)
preds <- data.table(x = seq(17, 90, 1), prediction = preds$fit, 
                    ymin = preds$fit - 1.96 * preds$se.fit,
                    ymax = preds$fit + 1.96 * preds$se.fit)

hosp %>%
  ggplot(aes(x = age_mid, y = hosp_rate, ymin = bottom, ymax = top)) +
  geom_ribbon(inherit.aes = FALSE, data = preds, aes(x = x, ymin = ymin, ymax = ymax), alpha = 0.4) +
  geom_line(inherit.aes = FALSE, data = preds, aes(x = x, y = prediction), col = "blue") +
  geom_point() +
  geom_errorbar() +
  labs(y = "Hospitalisation rate", x = "Age group mid-point") +
  cowplot::theme_cowplot()


hosp_rate <- data.table(age_grp = agelabs, mid = midpoints)
preds <- predict(hosp_fit, newdata = data.frame(age_mid = hosp_rate$mid), se.fit = TRUE)
hosp_rate$hosp_rate <- preds$fit / 100
hosp_rate$bottom_hosp_rate <- (preds$fit - 1.96 * preds$se.fit) / 100
hosp_rate$top_hosp_rate <- (preds$fit + 1.96 * preds$se.fit) / 100
hosp_rate[bottom_hosp_rate < 0, bottom_hosp_rate := 0]

hosp_rate <- merge(hosp_rate, hosp_age_prop, by = "age_grp")


# old
## Hospital acquired infection
# source("R/hospitalisation_rate.R")
# 
# hosp_acq_prop <- 0.11
# 
# hosp <- fr[type == "estimate" & variable == "infections" & location == "community"]
# 
# hosp[, location := "hospital"]
# hosp <- merge(hosp, hosp_rate, by = "age_grp")
# 
# cols <- c("bottom", "top", "lower", "upper", "median", "mean")
# 
# hosp[, (cols) := lapply(.SD, "*", hosp_acq_prop * hosp_rate * prop_of_hosp_acq), .SDcols = cols]
# 
# fr <- merge(fr, hosp_rate, by = "age_grp")
# fr[type == "estimate" & variable == "infections" & location == "community", (cols) := lapply(.SD, "*", 1 - (hosp_acq_prop * hosp_rate * prop_of_hosp_acq)), .SDcols = cols]
# 
# fr <- rbind(fr, hosp)
