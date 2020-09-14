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
# agebreaks <- c(0, 18, 25, 35, 45, 55, 65, 75, 100)
# agelabs <- c("0-17", "18-24", "25-34", "35-44", "45-54", "55-64", "65-74","75-100")