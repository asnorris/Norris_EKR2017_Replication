# locality-appendix.R: script for locality-level analyses presented in the
# appendix of the paper. If part of the script doesn't work, please first
# upgrade R and reinstall the packages loaded below. If the problem remains,
# please email ejdemyr@gmail.com. The script works as of R version 3.4.0
# (2017-04-21).

# Load packages
library(stargazer)
library(Hmisc)
library(sandwich)
library(tidyverse)
library(lfe)
library(broom)

# Load custom functions for analyses
source("functions.R")

e <- read_csv("locality-panel.csv")

# Add a dummy for second period and treated in second period; clean segregation
# dummies; add per capita measures; and add year for which public good was
# measured
e <- e %>%
  mutate(T = ifelse(period == "1999-2009", 1, 0),
         D = ifelse(T == 1 & mever == 1, 1, 0),
         seg3 = plyr::revalue(avediscat3, c(
           "[0.0425,0.401)" = "low",
           "[0.4010,0.490)" = "medium",
           "[0.4900,0.795]" = "high")),
         seg3 = factor(seg3, levels = c("low", "medium", "high")),
         bp1k_ln = log((nqborehole / total1k) + 0.001),
         qb_pc = nqborehole98 / total1k, #boreholes per capita in 1998
         year_pg = ifelse(period == "1999-2009", 2008, 1998)
         )

# Subsets
e_sub <- e %>%
  filter(match94 == 0,     #subset down to EAs that were not matched in 1994
         elf_ed >= 0.05)   #subset down to districts with baseline diversity

# Drop 160 EAs for which we don't have MP info
e_sub <- e_sub %>%
  filter(!D %in% NA) %>%
  group_by(eacode) %>%
  filter(n() == 2) %>% #this ensures that the 160 dropped observations are dropped for both time periods
  ungroup()

# Number of unique EAs
n_distinct(e_sub$eacode)


# Complete DiD results -------------

m1 <- felm(qb ~ T + D | eacode | 0 | eacode, data = e_sub, subset = seg3 == "low")
m2 <- felm(qb ~ T + D | eacode | 0 | eacode, data = e_sub, subset = seg3 == "medium")
m3 <- felm(qb ~ T + D | eacode | 0 | eacode, data = e_sub, subset = seg3 == "high")
m4 <- felm(qb ~ T + D + D:seg3 | eacode | 0 | eacode, data = e_sub)
m5 <- felm(qb ~ T + D + D:seg3 + qc + qs + pres_match | eacode | 0 | eacode, data = e_sub)
m6 <- felm(qb ~ T + D + D:seg3 + qc + qs + pres_match +
             T:urbancontains + T:log(popden+0.001) + T:elf +
             T:area_sqkm + T:distlilongwe_km + T:distblantyre_km +
             T:qb_pc | eacode | 0 | eacode, data = e_sub)
m7 <- felm(qb ~ T + D + D:avediskern + qc + qs + pres_match +
             T:urbancontains + T:log(popden+0.001) + T:elf +
             T:area_sqkm + T:distlilongwe_km + T:distblantyre_km +
             T:qb_pc | eacode | 0 | eacode, data = e_sub)

mods <- list(m1, m2, m3, m4, m5, m6, m7)

covlabs <- c(
  "Match with MP",
  "Match x Med. Seg.",
  "Match x High Seg.",
  "Match x Continuous Seg.",
  "Presence of Clinic",
  "Presence of School",
  "Presidential Ethnic Match",
  "Time",
  "Time x Urban",
  "Time x Pop. Density (log)",
  "Time x ELF",
  "Time x Land Area (sqkm)",
  "Time x Distance to Lilongwe",
  "Time x Distance to Blantyre",
  "Time x Boreholes per Capita, 1999"
)

edf <- e_sub %>% data.frame()
n_constit <- c(n_distinct(edf$constit[edf$seg3 == "low"]),
               n_distinct(edf$constit[edf$seg3 == "medium"]),
               n_distinct(edf$constit[edf$seg3 == "high"]),
               rep(n_distinct(edf$constit), 5)
               )

n_obs <- plyr::laply(mods, function(x) x$N)

newlines <- list(
  c("Locality Fixed Effects", rep("\\checkmark", 7)),
  c("Time Period Fixed Effects", rep("\\checkmark", 7)),
  c("Time Varying Controls", c("", "", "", "", "\\checkmark", "\\checkmark", "\\checkmark")),
  c("Fixed Controls x Time Period", c("", "", "", "", "", "\\checkmark", "\\checkmark")),
  c("Number of Electoral Districts", n_constit),
  c("Number of Localities", n_obs / 2),
  c("Number of Observations", n_obs)
)

collabs <- c("Low", "Med.", "High", "All", "All", "All", "All")

stargazer(mods,
          align = T,
          font.size = "small",
          title = "Segregation and Ethnic Favoritism in the Provision of Boreholes",
          label = "tab:ea_did_complete",
          covariate.labels = covlabs,
          order = c("D", "D:seg3medium", "D:seg3high", "D:avediskern", "qc", "qs", "pres_match"),
          column.labels = collabs,
          column.sep.width = "0pt",
          dep.var.labels = "Indicator for whether locality has a borehole",
          digits = 2,
          add.lines = newlines,
          omit.stat = c("adj.rsq", "rsq", "ser", "theta", "aic", "ll", "n"),
          notes.align = "l",
          notes.label = "",
          no.space = TRUE,
          out = "tables/ea_did_complete.tex")



# Including "all" districts (i.e., no ELF cutoff) -------------

e_all <- e %>%
  filter(match94 == 0) %>%
  filter(!D %in% NA) %>%
  group_by(eacode) %>%
  filter(n() == 2) %>%
  ungroup()

# Panel A
m1 <- felm(qb ~ T + D | eacode | 0 | eacode, data = e_all, subset = seg3 == "low")
m2 <- felm(qb ~ T + D | eacode | 0 | eacode, data = e_all, subset = seg3 == "medium")
m3 <- felm(qb ~ T + D | eacode | 0 | eacode, data = e_all, subset = seg3 == "high")
m4 <- felm(qb ~ T + D + D:seg3 | eacode | 0 | eacode, data = e_all)
m5 <- felm(qb ~ T + D + D:seg3 + qc + qs + pres_match | eacode | 0 | eacode, data = e_all)
m6 <- felm(qb ~ T + D + D:seg3 + qc + qs + pres_match +
             T:urbancontains + T:log(popden+0.001) + T:elf +
             T:area_sqkm + T:distlilongwe_km + T:distblantyre_km +
             T:qb_pc | eacode | 0 | eacode, data = e_all)
m7 <- felm(qb ~ T + D + D:avediskern + qc + qs + pres_match +
             T:urbancontains + T:log(popden+0.001) + T:elf +
             T:area_sqkm + T:distlilongwe_km + T:distblantyre_km +
             T:qb_pc | eacode | 0 | eacode, data = e_all)

mods <- list(m1, m2, m3, m4, m5, m6, m7)

covlabs <- c(
  "Match with MP",
  "Match x Med. Seg.",
  "Match x High Seg.",
  "Match x Continuous Seg."
)

edf <- e_all %>% data.frame()
n_constit <- c(n_distinct(edf$constit[edf$seg3 == "low"]),
               n_distinct(edf$constit[edf$seg3 == "medium"]),
               n_distinct(edf$constit[edf$seg3 == "high"]),
               rep(n_distinct(edf$constit), 5)
               )

n_obs <- plyr::laply(mods, function(x) x$N)

newlines <- list(
  c("Locality Fixed Effects", rep("\\checkmark", 7)),
  c("Time Period Fixed Effects", rep("\\checkmark", 7)),
  c("Time Varying Controls", c("", "", "", "", "\\checkmark", "\\checkmark", "\\checkmark")),
  c("Fixed Controls x Time Period", c("", "", "", "", "", "\\checkmark", "\\checkmark")),
  c("Number of Electoral Districts", n_constit),
  c("Number of Localities", n_obs / 2),
  c("Number of Observations", n_obs)
)

collabs <- c("Low", "Med.", "High", "All", "All", "All", "All")

stargazer(mods,
          omit = c("qc", "qs", "pres_match", "T", "T:urbancontains", "T:log(popden+0.001)" , "T:elf" , "T:area_sqkm" , "T:distlilongwe_km", "T:distblantyre_km", "T:qb_pc"),
          align = T,
          font.size = "small",
          title = "Segregation and Ethnic Favoritism in the Provision of Boreholes",
          label = "tab:did_all_a",
          covariate.labels = covlabs,
          column.labels = collabs,
          column.sep.width = "0pt",
          dep.var.labels = "Indicator for whether locality has a borehole",
          digits = 2,
          add.lines = newlines,
          omit.stat = c("adj.rsq", "rsq", "ser", "theta", "aic", "ll", "n"),
          notes.align = "l",
          notes.label = "",
          no.space = TRUE,
          out = "tables/ea_did_all_a.tex")


# Panel B
m1 <- felm(qb ~ T + matchshare | eacode | 0 | eacode, data = e_all, subset = seg3 == "low")
m2 <- felm(qb ~ T + matchshare | eacode | 0 | eacode, data = e_all, subset = seg3 == "medium")
m3 <- felm(qb ~ T + matchshare | eacode | 0 | eacode, data = e_all, subset = seg3 == "high")
m4 <- felm(qb ~ T + matchshare + matchshare:seg3 | eacode | 0 | eacode, data = e_all)
m5 <- felm(qb ~ T + matchshare + matchshare:seg3 + qc + qs + pres_match | eacode | 0 | eacode,
           data = e_all)
m6 <- felm(qb ~ T + matchshare + matchshare:seg3 + qc + qs + pres_match +
             T:urbancontains + T:log(popden+0.001) + T:elf + T:area_sqkm +
             T:distlilongwe_km + T:distblantyre_km + T:qb_pc | eacode | 0 |
             eacode, data = e_all)
m7 <- felm(qb ~ T + matchshare + matchshare:avediskern + qc + qs + pres_match +
             T:urbancontains + T:log(popden+0.001) + T:elf + T:area_sqkm +
             T:distlilongwe_km + T:distblantyre_km + T:qb_pc | eacode | 0 |
             eacode, data = e_all)

mods <- list(m1, m2, m3, m4, m5, m6, m7)

covlabs <- c(
  "Match with MP",
  "Match x Med. Seg.",
  "Match x High Seg.",
  "Match x Continuous Seg."
)

edf <- e_all %>% data.frame()
n_constit <- c(n_distinct(edf$constit[edf$seg3 == "low"]),
               n_distinct(edf$constit[edf$seg3 == "medium"]),
               n_distinct(edf$constit[edf$seg3 == "high"]),
               rep(n_distinct(edf$constit), 5)
               )

n_obs <- plyr::laply(mods, function(x) x$N)

newlines <- list(
  c("Locality Fixed Effects", rep("\\checkmark", 7)),
  c("Time Period Fixed Effects", rep("\\checkmark", 7)),
  c("Time Varying Controls", c("", "", "", "", "\\checkmark", "\\checkmark", "\\checkmark")),
  c("Fixed Controls x Time Period", c("", "", "", "", "", "\\checkmark", "\\checkmark")),
  c("Number of Electoral Districts", n_constit),
  c("Number of Localities", n_obs / 2),
  c("Number of Observations", n_obs)
)

collabs <- c("Low", "Med.", "High", "All", "All", "All", "All")

stargazer(mods,
          omit = c("qc", "qs", "pres_match", "T", "T:urbancontains", "T:log(popden+0.001)" , "T:elf" , "T:area_sqkm" , "T:distlilongwe_km", "T:distblantyre_km", "T:qb_pc"),
          align = T,
          font.size = "small",
          title = "Segregation and Ethnic Favoritism in the Provision of Boreholes",
          label = "tab:did_all_b",
          covariate.labels = covlabs,
          column.labels = collabs,
          column.sep.width = "0pt",
          dep.var.labels = "Indicator for whether locality has a borehole",
          digits = 2,
          add.lines = newlines,
          omit.stat = c("adj.rsq", "rsq", "ser", "theta", "aic", "ll", "n"),
          notes.align = "l",
          notes.label = "",
          no.space = TRUE,
          out = "tables/ea_did_all_b.tex")


# Including only rural electoral districts -------------

e_rural <- e %>%
  filter(match94 == 0,
         elf_ed >= 0.05,
         urbanprop_ed < 0.3,
         !D %in% NA) %>%
  group_by(eacode) %>%
  filter(n() == 2) %>%
  ungroup()

# Panel A
m1 <- felm(qb ~ T + D | eacode | 0 | eacode, data = e_rural, subset = seg3 == "low")
m2 <- felm(qb ~ T + D | eacode | 0 | eacode, data = e_rural, subset = seg3 == "medium")
m3 <- felm(qb ~ T + D | eacode | 0 | eacode, data = e_rural, subset = seg3 == "high")
m4 <- felm(qb ~ T + D + D:seg3 | eacode | 0 | eacode, data = e_rural)
m5 <- felm(qb ~ T + D + D:seg3 + qc + qs + pres_match | eacode | 0 | eacode, data = e_rural)
m6 <- felm(qb ~ T + D + D:seg3 + qc + qs + pres_match +
             T:urbancontains + T:log(popden+0.001) + T:elf +
             T:area_sqkm + T:distlilongwe_km + T:distblantyre_km +
             T:qb_pc | eacode | 0 | eacode, data = e_rural)
m7 <- felm(qb ~ T + D + D:avediskern + qc + qs + pres_match +
             T:urbancontains + T:log(popden+0.001) + T:elf +
             T:area_sqkm + T:distlilongwe_km + T:distblantyre_km +
             T:qb_pc | eacode | 0 | eacode, data = e_rural)

mods <- list(m1, m2, m3, m4, m5, m6, m7)

covlabs <- c(
  "Match with MP",
  "Match x Med. Seg.",
  "Match x High Seg.",
  "Match x Continuous Seg."
)

edf <- e_rural %>% data.frame()
n_constit <- c(n_distinct(edf$constit[edf$seg3 == "low"]),
               n_distinct(edf$constit[edf$seg3 == "medium"]),
               n_distinct(edf$constit[edf$seg3 == "high"]),
               rep(n_distinct(edf$constit), 5)
               )

n_obs <- plyr::laply(mods, function(x) x$N)

newlines <- list(
  c("Locality Fixed Effects", rep("\\checkmark", 7)),
  c("Time Period Fixed Effects", rep("\\checkmark", 7)),
  c("Time Varying Controls", c("", "", "", "", "\\checkmark", "\\checkmark", "\\checkmark")),
  c("Fixed Controls x Time Period", c("", "", "", "", "", "\\checkmark", "\\checkmark")),
  c("Number of Electoral Districts", n_constit),
  c("Number of Localities", n_obs / 2),
  c("Number of Observations", n_obs)
)

collabs <- c("Low", "Med.", "High", "All", "All", "All", "All")

stargazer(mods,
          omit = c("qc", "qs", "pres_match", "T", "T:urbancontains", "T:log(popden+0.001)" , "T:elf" , "T:area_sqkm" , "T:distlilongwe_km", "T:distblantyre_km", "T:qb_pc"),
          align = T,
          font.size = "small",
          title = "Segregation and Ethnic Favoritism in the Provision of Boreholes",
          label = "tab:did_rural_a",
          covariate.labels = covlabs,
          column.labels = collabs,
          column.sep.width = "0pt",
          dep.var.labels = "Indicator for whether locality has a borehole",
          digits = 2,
          add.lines = newlines,
          omit.stat = c("adj.rsq", "rsq", "ser", "theta", "aic", "ll", "n"),
          notes.align = "l",
          notes.label = "",
          no.space = TRUE,
          out = "tables/ea_did_rural_a.tex")


# Panel B
m1 <- felm(qb ~ T + matchshare | eacode | 0 | eacode, data = e_rural, subset = seg3 == "low")
m2 <- felm(qb ~ T + matchshare | eacode | 0 | eacode, data = e_rural, subset = seg3 == "medium")
m3 <- felm(qb ~ T + matchshare | eacode | 0 | eacode, data = e_rural, subset = seg3 == "high")
m4 <- felm(qb ~ T + matchshare + matchshare:seg3 | eacode | 0 | eacode, data = e_rural)
m5 <- felm(qb ~ T + matchshare + matchshare:seg3 + qc + qs + pres_match | eacode | 0 | eacode,
           data = e_rural)
m6 <- felm(qb ~ T + matchshare + matchshare:seg3 + qc + qs + pres_match +
             T:urbancontains + T:log(popden+0.001) + T:elf + T:area_sqkm +
             T:distlilongwe_km + T:distblantyre_km + T:qb_pc | eacode | 0 |
             eacode, data = e_rural)
m7 <- felm(qb ~ T + matchshare + matchshare:avediskern + qc + qs + pres_match +
             T:urbancontains + T:log(popden+0.001) + T:elf + T:area_sqkm +
             T:distlilongwe_km + T:distblantyre_km + T:qb_pc | eacode | 0 |
             eacode, data = e_rural)

mods <- list(m1, m2, m3, m4, m5, m6, m7)

covlabs <- c(
  "Match with MP",
  "Match x Med. Seg.",
  "Match x High Seg.",
  "Match x Continuous Seg."
)

edf <- e_rural %>% data.frame()
n_constit <- c(n_distinct(edf$constit[edf$seg3 == "low"]),
               n_distinct(edf$constit[edf$seg3 == "medium"]),
               n_distinct(edf$constit[edf$seg3 == "high"]),
               rep(n_distinct(edf$constit), 5)
               )

n_obs <- plyr::laply(mods, function(x) x$N)

newlines <- list(
  c("Locality Fixed Effects", rep("\\checkmark", 7)),
  c("Time Period Fixed Effects", rep("\\checkmark", 7)),
  c("Time Varying Controls", c("", "", "", "", "\\checkmark", "\\checkmark", "\\checkmark")),
  c("Fixed Controls x Time Period", c("", "", "", "", "", "\\checkmark", "\\checkmark")),
  c("Number of Electoral Districts", n_constit),
  c("Number of Localities", n_obs / 2),
  c("Number of Observations", n_obs)
)

collabs <- c("Low", "Med.", "High", "All", "All", "All", "All")

stargazer(mods,
          omit = c("qc", "qs", "pres_match", "T", "T:urbancontains", "T:log(popden+0.001)" , "T:elf" , "T:area_sqkm" , "T:distlilongwe_km", "T:distblantyre_km", "T:qb_pc"),
          align = T,
          font.size = "small",
          title = "Segregation and Ethnic Favoritism in the Provision of Boreholes",
          label = "tab:did_rural_b",
          covariate.labels = covlabs,
          column.labels = collabs,
          column.sep.width = "0pt",
          dep.var.labels = "Indicator for whether locality has a borehole",
          digits = 2,
          add.lines = newlines,
          omit.stat = c("adj.rsq", "rsq", "ser", "theta", "aic", "ll", "n"),
          notes.align = "l",
          notes.label = "",
          no.space = TRUE,
          out = "tables/ea_did_rural_b.tex")


# Excluding electoral districts in Machinga and Mangochi -------------

e_rm <- e %>%
  filter(match94 == 0,
         elf_ed >= 0.05,
         !admdistrict %in% c("Machinga", "Mangochi"),
         !D %in% NA) %>%
  group_by(eacode) %>%
  filter(n() == 2) %>%
  ungroup()

# Panel A
m1 <- felm(qb ~ T + D | eacode | 0 | eacode, data = e_rm, subset = seg3 == "low")
m2 <- felm(qb ~ T + D | eacode | 0 | eacode, data = e_rm, subset = seg3 == "medium")
m3 <- felm(qb ~ T + D | eacode | 0 | eacode, data = e_rm, subset = seg3 == "high")
m4 <- felm(qb ~ T + D + D:seg3 | eacode | 0 | eacode, data = e_rm)
m5 <- felm(qb ~ T + D + D:seg3 + qc + qs + pres_match | eacode | 0 | eacode, data = e_rm)
m6 <- felm(qb ~ T + D + D:seg3 + qc + qs + pres_match +
             T:urbancontains + T:log(popden+0.001) + T:elf +
             T:area_sqkm + T:distlilongwe_km + T:distblantyre_km +
             T:qb_pc | eacode | 0 | eacode, data = e_rm)
m7 <- felm(qb ~ T + D + D:avediskern + qc + qs + pres_match +
             T:urbancontains + T:log(popden+0.001) + T:elf +
             T:area_sqkm + T:distlilongwe_km + T:distblantyre_km +
             T:qb_pc | eacode | 0 | eacode, data = e_rm)

mods <- list(m1, m2, m3, m4, m5, m6, m7)

covlabs <- c(
  "Match with MP",
  "Match x Med. Seg.",
  "Match x High Seg.",
  "Match x Continuous Seg."
)

edf <- e_rm %>% data.frame()
n_constit <- c(n_distinct(edf$constit[edf$seg3 == "low"]),
               n_distinct(edf$constit[edf$seg3 == "medium"]),
               n_distinct(edf$constit[edf$seg3 == "high"]),
               rep(n_distinct(edf$constit), 5)
               )

n_obs <- plyr::laply(mods, function(x) x$N)

newlines <- list(
  c("Locality Fixed Effects", rep("\\checkmark", 7)),
  c("Time Period Fixed Effects", rep("\\checkmark", 7)),
  c("Time Varying Controls", c("", "", "", "", "\\checkmark", "\\checkmark", "\\checkmark")),
  c("Fixed Controls x Time Period", c("", "", "", "", "", "\\checkmark", "\\checkmark")),
  c("Number of Electoral Districts", n_constit),
  c("Number of Localities", n_obs / 2),
  c("Number of Observations", n_obs)
)

collabs <- c("Low", "Med.", "High", "All", "All", "All", "All")

stargazer(mods,
          omit = c("qc", "qs", "pres_match", "T", "T:urbancontains", "T:log(popden+0.001)" , "T:elf" , "T:area_sqkm" , "T:distlilongwe_km", "T:distblantyre_km", "T:qb_pc"),
          align = T,
          font.size = "small",
          title = "Segregation and Ethnic Favoritism in the Provision of Boreholes",
          label = "tab:did_rm_a",
          covariate.labels = covlabs,
          column.labels = collabs,
          column.sep.width = "0pt",
          dep.var.labels = "Indicator for whether locality has a borehole",
          digits = 2,
          add.lines = newlines,
          omit.stat = c("adj.rsq", "rsq", "ser", "theta", "aic", "ll", "n"),
          notes.align = "l",
          notes.label = "",
          no.space = TRUE,
          out = "tables/ea_did_rm_a.tex")


# Panel B
m1 <- felm(qb ~ T + matchshare | eacode | 0 | eacode, data = e_rm, subset = seg3 == "low")
m2 <- felm(qb ~ T + matchshare | eacode | 0 | eacode, data = e_rm, subset = seg3 == "medium")
m3 <- felm(qb ~ T + matchshare | eacode | 0 | eacode, data = e_rm, subset = seg3 == "high")
m4 <- felm(qb ~ T + matchshare + matchshare:seg3 | eacode | 0 | eacode, data = e_rm)
m5 <- felm(qb ~ T + matchshare + matchshare:seg3 + qc + qs + pres_match | eacode | 0 | eacode,
           data = e_rm)
m6 <- felm(qb ~ T + matchshare + matchshare:seg3 + qc + qs + pres_match +
             T:urbancontains + T:log(popden+0.001) + T:elf + T:area_sqkm +
             T:distlilongwe_km + T:distblantyre_km + T:qb_pc | eacode | 0 |
             eacode, data = e_rm)
m7 <- felm(qb ~ T + matchshare + matchshare:avediskern + qc + qs + pres_match +
             T:urbancontains + T:log(popden+0.001) + T:elf + T:area_sqkm +
             T:distlilongwe_km + T:distblantyre_km + T:qb_pc | eacode | 0 |
             eacode, data = e_rm)

mods <- list(m1, m2, m3, m4, m5, m6, m7)

covlabs <- c(
  "Match with MP",
  "Match x Med. Seg.",
  "Match x High Seg.",
  "Match x Continuous Seg."
)

edf <- e_rm %>% data.frame()
n_constit <- c(n_distinct(edf$constit[edf$seg3 == "low"]),
               n_distinct(edf$constit[edf$seg3 == "medium"]),
               n_distinct(edf$constit[edf$seg3 == "high"]),
               rep(n_distinct(edf$constit), 5)
               )

n_obs <- plyr::laply(mods, function(x) x$N)

newlines <- list(
  c("Locality Fixed Effects", rep("\\checkmark", 7)),
  c("Time Period Fixed Effects", rep("\\checkmark", 7)),
  c("Time Varying Controls", c("", "", "", "", "\\checkmark", "\\checkmark", "\\checkmark")),
  c("Fixed Controls x Time Period", c("", "", "", "", "", "\\checkmark", "\\checkmark")),
  c("Number of Electoral Districts", n_constit),
  c("Number of Localities", n_obs / 2),
  c("Number of Observations", n_obs)
)

collabs <- c("Low", "Med.", "High", "All", "All", "All", "All")

stargazer(mods,
          omit = c("qc", "qs", "pres_match", "T", "T:urbancontains", "T:log(popden+0.001)" , "T:elf" , "T:area_sqkm" , "T:distlilongwe_km", "T:distblantyre_km", "T:qb_pc"),
          align = T,
          font.size = "small",
          title = "Segregation and Ethnic Favoritism in the Provision of Boreholes",
          label = "tab:did_rm_b",
          covariate.labels = covlabs,
          column.labels = collabs,
          column.sep.width = "0pt",
          dep.var.labels = "Indicator for whether locality has a borehole",
          digits = 2,
          add.lines = newlines,
          omit.stat = c("adj.rsq", "rsq", "ser", "theta", "aic", "ll", "n"),
          notes.align = "l",
          notes.label = "",
          no.space = TRUE,
          out = "tables/ea_did_rm_b.tex")



# Robustness to different segregation cutoffs -------------

# Generate k different cutoffs
set.seed(22222)
k <- 15
segco <- data.frame(lseg = rep(NA, 15), useg = rep(NA, 15))
for(i in 1:k) {
  segco[i, 1] <- round(runif(1, 0.25, 0.4), 3)
  segco[i, 2] <- segco[i, 1] + 0.12
}
segco <- arrange(segco, lseg)

# Run DiD estimates for these different cutoffs
minseg <- min(e_sub$avediskern, na.rm = T)
maxseg <- max(e_sub$avediskern, na.rm = T)

estimates <- bind_rows(lapply(1:k, function(i) {

  df1 <- e_sub %>% filter(avediskern < segco[i, 1])
  df2 <- e_sub %>% filter(avediskern >= segco[i, 1], avediskern < segco[i, 2])
  df3 <- e_sub %>% filter(avediskern >= segco[i, 2])

  m1 <- felm(qb ~ T + D | eacode | 0 | eacode, data = df1)
  m2 <- felm(qb ~ T + D | eacode | 0 | eacode, data = df2)
  m3 <- felm(qb ~ T + D | eacode | 0 | eacode, data = df3)

  tidy_mods <- list(tidy(m1), tidy(m2), tidy(m3)) %>%
    bind_rows() %>%
    filter(term == "D") %>%
    mutate(ub = estimate + std.error * qnorm(0.95)) %>%
    mutate(lb = estimate - std.error * qnorm(0.95)) %>%
    mutate(segcat = c(paste0(minseg, "-", segco[i, 1]),
                      paste0(segco[i, 1], "-", segco[i, 2]),
                      paste0(segco[i, 2], "-", maxseg))) %>%
    mutate(seg = c("low", "medium", "high")) %>%
    mutate(seg = factor(seg, levels = c("low", "medium", "high"))) %>%
    mutate(group = i) %>%
    arrange(seg)

  return(tidy_mods)

}))

p <-
  ggplot(estimates, aes(x = seg, y = estimate)) +
  geom_errorbar(aes(ymin = lb, ymax = ub), alpha = 0.5, width = 0.15) +
  geom_point() +
  geom_line(aes(group = group), alpha = 0.5) +
  facet_wrap(~group, scale = "free", ncol=3) +
  geom_text(aes(label=segcat, y=ub+.02), size = 2.5) +
  geom_hline(yintercept = 0, lty=2) +
  xlab("\nSegregation Category") +
  ylab("DiD Estimate\n") +
  theme_classic() +
  theme(strip.background = element_blank()) +
  theme(strip.text.x = element_blank()) +
  theme(axis.text = element_text(size = 8)) +
  theme(axis.title = element_text(size = 11))


ggsave(
  "figures/did_borehole_cutoffs.pdf",
  plot = p,
  width = 8,
  height = 9.2,
  )


# Clustering at constituency-time period -------------

# Panel A
m1 <- felm(qb ~ T + D | eacode | 0 | constit+T, data = e_sub, subset = seg3 == "low")
m2 <- felm(qb ~ T + D | eacode | 0 | constit+T, data = e_sub, subset = seg3 == "medium")
m3 <- felm(qb ~ T + D | eacode | 0 | constit+T, data = e_sub, subset = seg3 == "high")
m4 <- felm(qb ~ T + D + D:seg3 | eacode | 0 | constit+T, data = e_sub)
m5 <- felm(qb ~ T + D + D:seg3 + qc + qs + pres_match | eacode | 0 | constit+T, data = e_sub)
m6 <- felm(qb ~ T + D + D:seg3 + qc + qs + pres_match  + T:urbancontains + T:log(popden+0.001) + T:elf + T:area_sqkm + T:distlilongwe_km + T:distblantyre_km + T:qb_pc | eacode | 0 |constit+T, data = e_sub)
m7 <- felm(qb ~ T + D + D:avediskern + qc + qs + pres_match +T:urbancontains + T:log(popden+0.001) + T:elf + T:area_sqkm + T:distlilongwe_km + T:distblantyre_km + T:qb_pc | eacode | 0 | constit+T, data = e_sub)

mods <- list(m1, m2, m3, m4, m5, m6, m7)

covlabs <- c(
    "Match with MP",
    "Match x Med. Seg.",
    "Match x High Seg.",
    "Match x Continuous Seg."
    )

edf <- e_sub %>% data.frame()
n_constit <- c(n_distinct(edf$constit[edf$seg3 == "low"]),
               n_distinct(edf$constit[edf$seg3 == "medium"]),
               n_distinct(edf$constit[edf$seg3 == "high"]),
               rep(n_distinct(edf$constit), 5)
               )

n_obs <- plyr::laply(mods, function(x) x$N)

newlines <- list(
    c("Locality Fixed Effects", rep("\\checkmark", 7)),
    c("Time Period Fixed Effects", rep("\\checkmark", 7)),
    c("Time Varying Controls", c("", "", "", "", "\\checkmark", "\\checkmark", "\\checkmark")),
    c("Fixed Controls x Time Period", c("", "", "", "", "", "\\checkmark", "\\checkmark")),
    c("Number of Electoral Districts", n_constit),
    c("Number of Localities", n_obs / 2),
    c("Number of Observations", n_obs)
    )

collabs <- c("Low", "Med.", "High", "All", "All", "All", "All")

stargazer(mods,
          omit = c("qc", "qs", "pres_match", "T",  "T:urbancontains", "T:log(popden+0.001)" , "T:elf" , "T:area_sqkm" , "T:distlilongwe_km", "T:distblantyre_km", "T:qb_pc"),
          align = T,
          font.size = "small",
          title = "Segregation and Ethnic Favoritism in the Provision of Boreholes (Standard Errors Clustered by Electoral District-Year)",
          label = "tab:did_se",
          covariate.labels = covlabs,
          column.labels = collabs,
          column.sep.width = "0pt",
          dep.var.labels = "Indicator for whether locality has a borehole",
          digits = 2,
          add.lines = newlines,
          omit.stat = c("adj.rsq", "rsq", "ser", "theta", "aic", "ll", "n"),
          notes.align = "l",
          notes.label = "",
          no.space = TRUE,
          out = "tables/ea_did_newSE.tex")

# Panel B
m1 <- felm(qb ~ T + matchshare | eacode | 0 | constit+T, data = e_sub, subset = seg3 == "low")
m2 <- felm(qb ~ T + matchshare | eacode | 0 | constit+T, data = e_sub, subset = seg3 == "medium")
m3 <- felm(qb ~ T + matchshare | eacode | 0 | constit+T, data = e_sub, subset = seg3 == "high")
m4 <- felm(qb ~ T + matchshare + matchshare:seg3 | eacode | 0 | constit+T, data = e_sub)
m5 <- felm(qb ~ T + matchshare + matchshare:seg3 + qc + qs + pres_match | eacode | 0 | constit+T, data = e_sub)
m6 <- felm(qb ~ T + matchshare + matchshare:seg3 + qc + qs + pres_match  + T:urbancontains + T:log(popden+0.001) + T:elf + T:area_sqkm + T:distlilongwe_km + T:distblantyre_km + T:qb_pc | eacode | 0 | constit+T, data = e_sub)
m7 <- felm(qb ~ T + matchshare + matchshare:avediskern + qc + qs + pres_match +T:urbancontains + T:log(popden+0.001) + T:elf + T:area_sqkm + T:distlilongwe_km + T:distblantyre_km + T:qb_pc | eacode | 0 | constit+T, data = e_sub)

mods <- list(m1, m2, m3, m4, m5, m6, m7)

stargazer(mods,
          omit = c("qc", "qs", "pres_match", "T",  "T:urbancontains", "T:log(popden+0.001)" , "T:elf" , "T:area_sqkm" , "T:distlilongwe_km", "T:distblantyre_km", "T:qb_pc"),
          align = T,
          font.size = "small",
          title = "Segregation and Ethnic Favoritism in the Provision of Boreholes (Standard Errors Clustered by Electoral District-Year)",
          label = "tab:did_contiv_newSE",
          covariate.labels = covlabs,
          column.labels = collabs,
          column.sep.width = "0pt",
          dep.var.labels = "Indicator for whether locality has a borehole",
          digits = 2,
          add.lines = newlines,
          omit.stat = c("adj.rsq", "rsq", "ser", "theta", "aic", "ll", "n"),
          notes.align = "l",
          notes.label = "",
          no.space = TRUE,
          out = "tables/ea_did_contiv_newSE.tex")


# Schools -------------

# Binary match
m1 <- felm(qs ~ T + D | eacode | 0 | eacode, data = e_sub, subset = seg3 == "low")
m2 <- felm(qs ~ T + D | eacode | 0 | eacode, data = e_sub, subset = seg3 == "medium")
m3 <- felm(qs ~ T + D | eacode | 0 | eacode, data = e_sub, subset = seg3 == "high")
m4 <- felm(qs ~ T + D + D:seg3 | eacode | 0 | eacode, data = e_sub)
m5 <- felm(qs ~ T + D + D:seg3 + qb + qc + pres_match | eacode | 0 | eacode, data = e_sub)
m6 <- felm(qs ~ T + D + D:avediskern + qb + qc + pres_match | eacode | 0 | eacode, data = e_sub)

mods <- list(m1, m2, m3, m4, m5, m6)

covlabs <- c(
    "Ethnic Match with MP",
    "Ethnic Match x Med. Segregation",
    "Ethnic Match x High Segregation",
    "Ethnic Match x Continuous Segregation"
    )

edf <- e_sub %>% data.frame()
n_constit <- c(n_distinct(edf$constit[edf$seg3 == "low"]),
               n_distinct(edf$constit[edf$seg3 == "medium"]),
               n_distinct(edf$constit[edf$seg3 == "high"]),
               rep(n_distinct(edf$constit), 3)
               )

n_eas <- plyr::laply(mods, function(x) x$N)

newlines <- list(
    c("Locality Fixed Effects", rep("\\checkmark", 6)),
    c("Time Period Fixed Effects", rep("\\checkmark", 6)),
    c("Time Varying Controls", c("", "", "", "", "\\checkmark", "\\checkmark")),
    c("Number of Constituencies", n_constit),
    c("Number of Localities", n_eas)
    )

collabs <- c("Low", "Med.", "High", "All", "All", "All")

stargazer(mods,
          omit = c("qb", "qc", "pres_match", "T"),
          align = T,
          font.size = "small",
          title = "Segregation and Ethnic Favoritism in the Provision of Schools (Binary Match)",
          label = "tab:ea_did_sch",
          covariate.labels = covlabs,
          column.labels = collabs,
          column.sep.width = "0pt",
          dep.var.labels = "Indicator for presence of a school",
          digits = 2,
          add.lines = newlines,
          omit.stat = c("rsq", "ser", "theta", "aic", "ll", "n"),
          notes.align = "l",
          notes.label = "",
          no.space = TRUE,
          out = "tables/ea_did_sch.tex")


# Proportion coethnic
m1 <- felm(qs ~ T + matchshare | eacode | 0 | eacode, data = e_sub, subset = seg3 == "low")
m2 <- felm(qs ~ T + matchshare | eacode | 0 | eacode, data = e_sub, subset = seg3 == "medium")
m3 <- felm(qs ~ T + matchshare | eacode | 0 | eacode, data = e_sub, subset = seg3 == "high")
m4 <- felm(qs ~ T + matchshare + matchshare:seg3 | eacode | 0 | eacode, data = e_sub)
m5 <- felm(qs ~ T + matchshare + matchshare:seg3 + qb + qc + pres_match | eacode | 0 | eacode, data = e_sub)
m6 <- felm(qs ~ T + matchshare + matchshare:avediskern + qb + qc + pres_match | eacode | 0 | eacode, data = e_sub)

mods <- list(m1, m2, m3, m4, m5, m6)

covlabs <- c(
    "Ethnic Match with MP",
    "Ethnic Match x Med. Segregation",
    "Ethnic Match x High Segregation",
    "Ethnic Match x Continuous Segregation"
    )

edf <- e_sub %>% data.frame()
n_constit <- c(n_distinct(edf$constit[edf$seg3 == "low"]),
               n_distinct(edf$constit[edf$seg3 == "medium"]),
               n_distinct(edf$constit[edf$seg3 == "high"]),
               rep(n_distinct(edf$constit), 3)
               )

n_eas <- plyr::laply(mods, function(x) x$N)

newlines <- list(
    c("Locality Fixed Effects", rep("\\checkmark", 6)),
    c("Time Period Fixed Effects", rep("\\checkmark", 6)),
    c("Time Varying Controls", c("", "", "", "", "\\checkmark", "\\checkmark")),
    c("Number of Constituencies", n_constit),
    c("Number of Localities", n_eas)
    )

collabs <- c("Low", "Med.", "High", "All", "All", "All")

stargazer(mods,
          omit = c("qb", "qc", "pres_match", "T"),
          align = T,
          font.size = "small",
          title = "Segregation and Ethnic Favoritism in the Provision of Schools (Proportion Coethnic)",
          label = "tab:ea_did_sch_contiv",
          covariate.labels = covlabs,
          column.labels = collabs,
          column.sep.width = "0pt",
          dep.var.labels = "Indicator for presence of a school",
          digits = 2,
          add.lines = newlines,
          omit.stat = c("rsq", "ser", "theta", "aic", "ll", "n"),
          notes.align = "l",
          notes.label = "",
          no.space = TRUE,
          out = "tables/ea_did_sch_contiv.tex")



# Clinics -------------

# Binary match
m1 <- felm(qc ~ T + D | eacode | 0 | eacode, data = e_sub, subset = seg3 == "low")
m2 <- felm(qc ~ T + D | eacode | 0 | eacode, data = e_sub, subset = seg3 == "medium")
m3 <- felm(qc ~ T + D | eacode | 0 | eacode, data = e_sub, subset = seg3 == "high")
m4 <- felm(qc ~ T + D + D:seg3 | eacode | 0 | eacode, data = e_sub)
m5 <- felm(qc ~ T + D + D:seg3 + qb + qs + pres_match | eacode | 0 | eacode, data = e_sub)
m6 <- felm(qc ~ T + D + D:avediskern + qb + qs + pres_match | eacode | 0 | eacode, data = e_sub)

mods <- list(m1, m2, m3, m4, m5, m6)

covlabs <- c(
    "Ethnic Match with MP",
    "Ethnic Match x Med. Segregation",
    "Ethnic Match x High Segregation",
    "Ethnic Match x Continuous Segregation"
    )

edf <- e_sub %>% data.frame()
n_constit <- c(n_distinct(edf$constit[edf$seg3 == "low"]),
               n_distinct(edf$constit[edf$seg3 == "medium"]),
               n_distinct(edf$constit[edf$seg3 == "high"]),
               rep(n_distinct(edf$constit), 3)
               )

n_eas <- plyr::laply(mods, function(x) x$N)

newlines <- list(
    c("Locality Fixed Effects", rep("\\checkmark", 6)),
    c("Time Period Fixed Effects", rep("\\checkmark", 6)),
    c("Time Varying Controls", c("", "", "", "", "\\checkmark", "\\checkmark")),
    c("Number of Constituencies", n_constit),
    c("Number of Localities", n_eas)
    )

collabs <- c("Low", "Med.", "High", "All", "All", "All")

stargazer(mods,
          omit = c("qb", "qs", "pres_match", "T"),
          align = T,
          font.size = "small",
          title = "Segregation and Ethnic Favoritism in the Provision of Clinics (Binary Match)",
          label = "tab:ea_did_cl",
          covariate.labels = covlabs,
          column.labels = collabs,
          column.sep.width = "0pt",
          dep.var.labels = "Indicator for presence of a clinic",
          digits = 2,
          add.lines = newlines,
          omit.stat = c("rsq", "ser", "theta", "aic", "ll", "n"),
          notes.align = "l",
          notes.label = "",
          no.space = TRUE,
          out = "tables/ea_did_cl.tex")


# Proportion coethnic
m1 <- felm(qc ~ T + matchshare | eacode | 0 | eacode, data = e_sub, subset = seg3 == "low")
m2 <- felm(qc ~ T + matchshare | eacode | 0 | eacode, data = e_sub, subset = seg3 == "medium")
m3 <- felm(qc ~ T + matchshare | eacode | 0 | eacode, data = e_sub, subset = seg3 == "high")
m4 <- felm(qc ~ T + matchshare + matchshare:seg3 | eacode | 0 | eacode, data = e_sub)
m5 <- felm(qc ~ T + matchshare + matchshare:seg3 + qb + qs + pres_match | eacode | 0 | eacode, data = e_sub)
m6 <- felm(qc ~ T + matchshare + matchshare:avediskern + qb + qs + pres_match | eacode | 0 | eacode, data = e_sub)

mods <- list(m1, m2, m3, m4, m5, m6)

covlabs <- c(
    "Ethnic Match with MP",
    "Ethnic Match x Med. Segregation",
    "Ethnic Match x High Segregation",
    "Ethnic Match x Continuous Segregation"
    )

edf <- e_sub %>% data.frame()
n_constit <- c(n_distinct(edf$constit[edf$seg3 == "low"]),
               n_distinct(edf$constit[edf$seg3 == "medium"]),
               n_distinct(edf$constit[edf$seg3 == "high"]),
               rep(n_distinct(edf$constit), 3)
               )

n_eas <- plyr::laply(mods, function(x) x$N)

newlines <- list(
    c("Locality Fixed Effects", rep("\\checkmark", 6)),
    c("Time Period Fixed Effects", rep("\\checkmark", 6)),
    c("Time Varying Controls", c("", "", "", "", "\\checkmark", "\\checkmark")),
    c("Number of Constituencies", n_constit),
    c("Number of Localities", n_eas)
    )

collabs <- c("Low", "Med.", "High", "All", "All", "All")

stargazer(mods,
          omit = c("qb", "qs", "pres_match", "T"),
          align = T,
          font.size = "small",
          title = "Segregation and Ethnic Favoritism in the Provision of Clinics (Continuous Match)",
          label = "tab:ea_did_cl_contiv",
          covariate.labels = covlabs,
          column.labels = collabs,
          column.sep.width = "0pt",
          dep.var.labels = "Indicator for presence of a clinic",
          digits = 2,
          add.lines = newlines,
          omit.stat = c("rsq", "ser", "theta", "aic", "ll", "n"),
          notes.align = "l",
          notes.label = "",
          no.space = TRUE,
          out = "tables/ea_did_cl_contiv.tex")


# Cross-sectional EA analysis -------------
ea <- read.csv("locality.csv")

# Ensure constituency is factor
ea <- ea %>% mutate(constit = factor(constit))

ea$lg_area <- log(ea$area_sqkm+0.001)
ea$lg_popden <- log(ea$popden+0.001)

# Set reference category for fixed effect to Machinga North East (~median intercept)
ea$constit <- relevel(ea$constit, ref = "Machinga North East")

# Boreholes per capita
ea$nqb_pc98 <- ea$nqborehole98/ea$total1k

f1 <- Formula(qnewborehole ~ mever + mever:avediskern + elf + match_pres + total1k + lg_popden + nqb_pc98 + urbancontains  + distlilongwe_km + distblantyre_km + lg_area + constit)
f2 <- Formula(qnewborehole ~ mever + mever:as.factor(avediscat3) + elf + match_pres + total1k + lg_popden+ nqb_pc98 + urbancontains  + distlilongwe_km + distblantyre_km + lg_area + constit)
f3 <- Formula(qnewborehole ~ avematchshare + avematchshare:avediskern + elf + match_pres + total1k + lg_popden + nqb_pc98 + urbancontains  + distlilongwe_km + distblantyre_km + lg_area + constit)
f4 <- Formula(qnewborehole ~  avematchshare + avematchshare:as.factor(avediscat3)+ elf + match_pres + total1k + lg_popden + nqb_pc98 + urbancontains  + distlilongwe_km + distblantyre_km + lg_area + constit)

m1 <- glm(f1, family=binomial(link="logit"), data=ea[ea$elf_ed > 0.05,])
m2 <- glm(f2, family=binomial(link="logit"), data=ea[ea$elf_ed > 0.05,])
m3 <- glm(f3, family=binomial(link="logit"), data=ea[ea$elf_ed > 0.05,])
m4 <- glm(f4, family=binomial(link="logit"), data=ea[ea$elf_ed > 0.05,])

covlabs <- c(
  "Ethnic Match with MP",
  "Percent Coethnic with MP",
  "Ethnic Match x Continuous Segregation",
  "Ethnic Match x Med. Segregation",
  "Ethnic Match x High Segregation",
  "Percent Coethnic with MP x Continuous Segregation",
  "Percent Coethnic with MP x Med. Segregation",
  "Percent Coethnic with MP x High Segregation"
)

newlines <- list(
  c("Constituency Fixed Effects", rep("\\checkmark", 4)),
  c("Control Variables", rep("\\checkmark", 4))
)

stargazer(m1, m2, m3, m4,
          omit = c("constit", "total1k" , "lg_popden" , "nqb_pc98" , "urbancontains" , "distlilongwe_km",  "distblantyre_km", "elf", "lg_area",  "match_pres"),
          align = T,
          font.size = "small",
          add.lines = newlines,
          title = "Segregation and Ethnic Favoritism in the Provision of Boreholes",
          label = "tab:ea_cross",
          covariate.labels = covlabs,
          column.sep.width = "0pt",
          dep.var.labels = "New Borehole",
          digits = 2,
          omit.stat = c("rsq", "ser", "theta", "aic", "ll", "n"),
          notes.align = "l",
          notes.label = "",
          no.space = TRUE,
          out = "tables/ea_cross.tex")
