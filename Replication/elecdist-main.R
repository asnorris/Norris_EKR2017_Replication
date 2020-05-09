# elecdist-main.R: script for election district analyses presented in the body
# of the paper. If part of the script doesn't work, please first upgrade R and
# reinstall the packages loaded below. If the problem remains, please email
# ejdemyr@gmail.com. The script works as of R version 3.4.0 (2017-04-21).

# Load packages
## Comment out package 'Hmisc' because it won't run
library(stargazer)
# library(Hmisc)
library(sandwich)
library(tidyverse)

# Load custom functions for analyses
source("Replication/functions.R")

# Read in data
ed <- read.csv("Replication/elecdist.csv")

# Ensure admin. district variable is factor
ed <- ed %>% mutate(admdistrict = factor(admdistrict))

# Set reference category for fixed effect to Chikwawa (~median intercept)
ed$admdistrict <- relevel(ed$admdistrict, ref = "Chikwawa")

## Create Table 1. Table 1 is the regression output of regressing a number of
## variables including segregation against the change in the number of boreholes
## between 1998 and 2008 

# Table 1 -------------
f1 <- formula(cqborehole ~ seg + elf_ed + I(log(popden)) + I(log(urbanprop + 0.001)) +
                log(area_sqkm) + nqborehole98_p10k + admdistrict)

f2 <- formula(cqborehole ~ seg + elf_ed + I(log(popden)) + I(log(urbanprop + 0.001)) +
                log(area_sqkm) + nqborehole98_p10k + avgwinmargpct + coethnic_share +
                prescoethavg + I(log(dist_closest + 1)) + ngowatern_p10k +
                ngoalln_p10k + admdistrict)

# Use function ed_mod to run the formula above for continuous
# segregation and dummied segregation
mods1 <- lapply(c("avediskern", "avediscat3"), function(seg) {
  ed_mod(df = ed, formula = f1, seg = seg, elfco = 0.05)
})

mods2 <- lapply(c("avediskern", "avediscat3"), function(seg) {
  ed_mod(df = ed, formula = f2, seg = seg, elfco = 0.05)
})

covlabs <- c("Segregation (continuous)",
             "Dummy for Medium Segregation",
             "Dummy for High Segregation",
             "Ethnic Diversity (ELF)",
             "Population Density (ln)",
             "Urban Proportion (ln)",
             "Land Area (square KM) (ln)",
             "Boreholes per 10,000 residents in 1998",
             "Electoral Competitiveness",
             "MP Coethnic Population Share",
             "President Coethnic Population Share",
             "Distance to Nearest City (ln)",
             "Water Aid Projects per 10,000 residents",
             "All Aid Projects per 10,000 residents"
             )

stargazer(mods1[1], mods2[1], mods1[2], mods2[2],
          omit = c("admdistrict"),
          align = T,
          font.size = "small",
          title = "Segregation and Borehole Investments across Electoral Districts",
          label = "tab:ed_boreholes",
          covariate.labels = covlabs,
          dep.var.labels = "Number of New Boreholes",
          digits = 2,
          add.lines = list(c("Admin. District Fixed Effects", rep("\\checkmark", 4))),
          omit.stat = c("rsq", "ser", "theta", "aic", "ll"),
          notes.align = "l",
          notes.label = "",
          no.space = TRUE,
          out = "Replication/tables/ed_boreholes.tex")


## Simulate the data to find confidence intervals

# Simulations of expected values and confidence intervals -------------
set.seed(222)
m <- 10000

# Continuous segregation (model 1)
seg_vals <- quantile(ed$avediskern, probs = c(0.1, 0.9), na.rm = T)
n_admdist <- sum(grepl("admdist", names(coef(mods1[[1]]))))

xset_low <- c("(Intercept)" = 1,
              "seg" = seg_vals[1],
              "elf_ed" = mean(ed$elf_ed),
              "popden" = mean(log(ed$popden)),
              "urbanprop" = mean(log(ed$urbanprop + 0.001)),
              "area_sqkm" = mean(log(ed$area_sqkm)),
              "nqborehole98_p10k" = mean(ed$nqborehole98_p10k),
              rep(0, n_admdist))

xset_high <- c("(Intercept)" = 1,
               "seg" = seg_vals[2],
               "elf_ed" = mean(ed$elf_ed),
               "popden" = mean(log(ed$popden)),
               "urbanprop" = mean(log(ed$urbanprop + 0.001)),
               "area_sqkm" = mean(log(ed$area_sqkm)),
               "nqborehole98_p10k" = mean(ed$nqborehole98_p10k),
               rep(0, n_admdist))

vcov_b <- ed_se(mod = mods1[[1]], output = "regvcov")
ed_sim(m = m, model = mods1[[1]], vcov = vcov_b, xset_low, xset_high)

# Compare with ELF (model 1)
elf_vals <- quantile(ed$elf_ed, probs = c(0.1, 0.9), na.rm = T)

xset_low <- c("(Intercept)" = 1,
              "seg" = mean(ed$avediskern),
              "elf_ed" = elf_vals[1],
              "popden" = mean(log(ed$popden)),
              "urbanprop" = mean(log(ed$urbanprop + 0.001)),
              "area_sqkm" = mean(log(ed$area_sqkm)),
              "nqborehole98_p10k" = mean(ed$nqborehole98_p10k),
              rep(0, n_admdist))

xset_high <- c("(Intercept)" = 1,
               "seg" = mean(ed$avediskern),
               "elf_ed" = elf_vals[2],
               "popden" = mean(log(ed$popden)),
               "urbanprop" = mean(log(ed$urbanprop + 0.001)),
               "area_sqkm" = mean(log(ed$area_sqkm)),
               "nqborehole98_p10k" = mean(ed$nqborehole98_p10k),
               rep(0, n_admdist))

ed_sim(m = m, model = mods1[[1]], vcov = vcov_b, xset_low, xset_high)

# Dummy variables for segregation (model 3)
xset_low <- c("(Intercept)" = 1,
              "seg_med" = 0,
              "seg_high" = 0,
              "elf_ed" = mean(ed$elf_ed),
              "popden" = mean(log(ed$popden)),
              "urbanprop" = mean(log(ed$urbanprop + 0.001)),
              "area_sqkm" = mean(log(ed$area_sqkm)),
              "nqborehole98_p10k" = mean(ed$nqborehole98_p10k),
              rep(0, n_admdist))

xset_med <- c("(Intercept)" = 1,
              "seg_med" = 1,
              "seg_high" = 0,
              "elf_ed" = mean(ed$elf_ed),
              "popden" = mean(log(ed$popden)),
              "urbanprop" = mean(log(ed$urbanprop + 0.001)),
              "area_sqkm" = mean(log(ed$area_sqkm)),
              "nqborehole98_p10k" = mean(ed$nqborehole98_p10k),
              rep(0, n_admdist))

xset_high <- c("(Intercept)" = 1,
               "seg_med" = 0,
               "seg_high" = 1,
               "elf_ed" = mean(ed$elf_ed),
               "popden" = mean(log(ed$popden)),
               "urbanprop" = mean(log(ed$urbanprop + 0.001)),
               "area_sqkm" = mean(log(ed$area_sqkm)),
               "nqborehole98_p10k" = mean(ed$nqborehole98_p10k),
               rep(0, n_admdist))

vcov_b2 <- ed_se(mod = mods1[[2]], output = "regvcov")

ed_sim(m = m, model = mods1[[2]], vcov = vcov_b2, xset_low, xset_med)
ed_sim(m = m, model = mods1[[2]], vcov = vcov_b2, xset_low, xset_high)
