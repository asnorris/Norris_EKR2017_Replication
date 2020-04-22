# locality-main.R: script for locality-level analyses presented in the body of
# the paper. If part of the script doesn't work, please first upgrade R and
# reinstall the packages loaded below. If the problem remains, please email
# ejdemyr@gmail.com. The script works as of R version 3.4.0 (2017-04-21).

# Load packages
library(stargazer)
library(Hmisc)
library(sandwich)
library(tidyverse)
library(lfe)

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

# Check change in number of boreholes
e_change <- e_sub %>%
  arrange(eacode, period) %>%
  group_by(eacode) %>%
  summarise(ch = nqborehole[2] - nqborehole[1])

prop.table(table(e_change$ch > 1))
with(e_change %>% filter(ch > 1), prop.table(table(ch)))

# What % had at least one borehole in 1998?
with(e_sub %>% filter(period == "1994-1999"), prop.table(table(qb == 1)))


# Figure 3 -------------
plot_dta <- e_sub %>%
  group_by(seg3, year_pg, mever) %>%
  summarise(y = mean(qb)) %>%
  spread(year_pg, y) %>%
  data.frame()

plot_dta

x <- c(1, 2)
plot_dta <- plot_dta[, 3:4]

pdf("figures/did_boreholes.pdf", width=7, height=4)
par(mfrow=c(1,3))
plot(x, plot_dta[1,],
     pch=19, xlim=c(.5,2.5), ylim=c(0, 0.55), las=1, cex=1.5, xaxt="n",
        main="Low Segregation", ylab="Pr(Locality Has Borehole)", xlab="")
	points(x, plot_dta[2, ], pch=24, cex=1.5)
	segments(x[1], plot_dta[1,1], x[2], plot_dta[1,2], lty=2)
	segments(x[1], plot_dta[2,1], x[2], plot_dta[2,2], lty=2)
	axis(1, at=x, labels=c("1998", "2008"))
	legend(.5, 55/100, legend=c("Coethnic after 1998", "Never Coethnic"),
               bty = "n", pch=c(24, 19), cex=0.85)

plot(x, plot_dta[3,], pch=19, xlim=c(.5,2.5), ylim=c(0, 0.55), las=1, cex=1.5, xaxt="n",
        main="Medium Segregation", ylab="", xlab="")
	points(x, plot_dta[4, ], pch=24, cex=1.5)
	segments(x[1], plot_dta[3,1], x[2], plot_dta[3,2], lty=2)
	segments(x[1], plot_dta[4,1], x[2], plot_dta[4,2], lty=2)
	axis(1, at=x, labels=c("1998", "2008"))


plot(x, plot_dta[5,], pch=19, xlim=c(.5,2.5), ylim=c(0, 0.55), las=1, cex=1.5, xaxt="n",
        main="High Segregation", ylab="", xlab="")
	points(x, plot_dta[6, ], pch=24, cex=1.5)
	segments(x[1], plot_dta[5,1], x[2], plot_dta[5,2], lty=2)
	segments(x[1], plot_dta[6,1], x[2], plot_dta[6,2], lty=2)
	axis(1, at=x, labels=c("1998", "2008"))
dev.off()


# Table 2, Panel A -------------

# Across three levels of segregation
m1 <- felm(qb ~ T + D | eacode | 0 | eacode, data = e_sub, subset = seg3 == "low")
m2 <- felm(qb ~ T + D | eacode | 0 | eacode, data = e_sub, subset = seg3 == "medium")
m3 <- felm(qb ~ T + D | eacode | 0 | eacode, data = e_sub, subset = seg3 == "high")

# Interact instead
m4 <- felm(qb ~ T + D + D:seg3 | eacode | 0 | eacode, data = e_sub)

# Time-varying controls
m5 <- felm(qb ~ T + D + D:seg3 + qc + qs + pres_match | eacode | 0 | eacode, data = e_sub)

# Time-varying controls + invariant controls interacted with time dummies
m6 <- felm(qb ~ T + D + D:seg3 + qc + qs + pres_match +
             T:urbancontains + T:log(popden+0.001) + T:elf +
             T:area_sqkm + T:distlilongwe_km + T:distblantyre_km +
             T:qb_pc | eacode | 0 | eacode, data = e_sub)

# Continuous segregation: Time-varying controls + invariant controls interacted with time dummies
m7 <- felm(qb ~ T + D + D:avediskern + qc + qs + pres_match +
             T:urbancontains + T:log(popden+0.001) + T:elf +
             T:area_sqkm + T:distlilongwe_km + T:distblantyre_km +
             T:qb_pc | eacode | 0 | eacode, data = e_sub)

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
          omit = c("qc", "qs", "pres_match", "T", "T:urbancontains", "T:log(popden+0.001)" , "T:elf" , "T:area_sqkm" , "T:distlilongwe_km", "T:distblantyre_km", "T:qb_pc"),
          align = T,
          font.size = "small",
          title = "Segregation and Ethnic Favoritism in the Provision of Boreholes",
          label = "tab:did",
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
          out = "tables/ea_did_a.tex")


# Table 2, Panel B -------------
m1 <- felm(qb ~ T + matchshare | eacode | 0 | eacode, data = e_sub, subset = seg3 == "low")
m2 <- felm(qb ~ T + matchshare | eacode | 0 | eacode, data = e_sub, subset = seg3 == "medium")
m3 <- felm(qb ~ T + matchshare | eacode | 0 | eacode, data = e_sub, subset = seg3 == "high")

# Interact instead
m4 <- felm(qb ~ T + matchshare + matchshare:seg3 | eacode | 0 | eacode, data = e_sub)

# Time-varying controls
m5 <- felm(qb ~ T + matchshare + matchshare:seg3 + qc + qs + pres_match | eacode | 0 | eacode,
           data = e_sub)

# Time-varying controls + invariant controls interacted with time dummies
m6 <- felm(qb ~ T + matchshare + matchshare:seg3 + qc + qs + pres_match +
             T:urbancontains + T:log(popden+0.001) + T:elf + T:area_sqkm +
             T:distlilongwe_km + T:distblantyre_km + T:qb_pc | eacode | 0 |
             eacode, data = e_sub)

# Time-varying controls + invariant controls interacted with time dummieswith continuous segregation
m7 <- felm(qb ~ T + matchshare + matchshare:avediskern + qc + qs + pres_match +
             T:urbancontains + T:log(popden+0.001) + T:elf + T:area_sqkm +
             T:distlilongwe_km + T:distblantyre_km + T:qb_pc | eacode | 0 |
             eacode, data = e_sub)

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
          omit = c("qc", "qs", "pres_match", "T", "T:urbancontains", "T:log(popden+0.001)" , "T:elf" , "T:area_sqkm" , "T:distlilongwe_km", "T:distblantyre_km", "T:qb_pc"),
          align = T,
          font.size = "small",
          title = "Segregation and Ethnic Favoritism in the Provision of Boreholes",
          label = "tab:did_contiv",
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
          out = "tables/ea_did_b.tex")
