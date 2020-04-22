# elecdist-appendix.R: script for election district analyses presented in the
# appendix of the paper. If part of the script doesn't work, please first
# upgrade R and reinstall the packages loaded below. If the problem remains,
# please email ejdemyr@gmail.com. The script works as of R version 3.4.0
# (2017-04-21).

# Load packages
library(stargazer)
# library(Hmisc)
library(sandwich)
library(tidyverse)

## edited file path to allow function to work

# Load custom functions for analyses
source("Replication/functions.R")

## edited file path to allow function to work

# Read in data
ed <- read.csv("Replication/elecdist.csv")

# Ensure admin. district variable is factor
ed <- ed %>% mutate(admdistrict = factor(admdistrict))

# Set reference category for fixed effect to Chikwawa (~median intercept)
ed$admdistrict <- relevel(ed$admdistrict, ref = "Chikwawa")

# Scatter plot: segregation v new boreholes -------------
plot_df1 <- ed %>%
  filter(elf_ed >= 0.05) %>%
  mutate(pop1k = pop / 1000)

plot_df2 <- plot_df1 %>%
  filter(cqborehole <= mean(cqborehole) + 2 * sd(cqborehole))

# Plot
p <- ggplot(plot_df1, aes(x = avediskern, y = cqborehole)) +
  geom_point(aes(size = pop1k)) +
  geom_smooth(aes(weight = pop1k), method = "loess", se = F) +
  geom_smooth(data = plot_df2, lty = 2,
              aes(x = avediskern, y = cqborehole, weight = pop1k),
              se = F, method = "loess") +
  scale_size("Population\n(1,000s)", range = c(0, 5)) +
  theme_classic() +
  xlab("\nSegregation") +
  ylab("Number of New Boreholes\n")  +
  theme_smpl

ggsave(
  "figures/scatter_constit.pdf",
  plot = p,
  width = 7.5,
  height = 4,
  )


# Including all electoral districts -------------

f1 <- formula(cqborehole ~ seg + elf_ed + I(log(popden)) + I(log(urbanprop + 0.001)) +
                log(area_sqkm) + nqborehole98_p10k + admdistrict)

f2 <- formula(cqborehole ~ seg + elf_ed + I(log(popden)) + I(log(urbanprop + 0.001)) +
                log(area_sqkm) + nqborehole98_p10k + avgwinmargpct + coethnic_share +
                prescoethavg + I(log(dist_closest + 1)) + ngowatern_p10k +
                ngoalln_p10k + admdistrict)

# Use function ed_mod to run the formula above for continuous
# segregation and dummied segregation
mods1 <- lapply(c("avediskern", "avediscat3"), function(seg) {
  ed_mod(df = ed, formula = f1, seg = seg, elfco = 0)
})

mods2 <- lapply(c("avediskern", "avediscat3"), function(seg) {
  ed_mod(df = ed, formula = f2, seg = seg, elfco = 0)
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
          title = "Segregation and Borehole Investments across Electoral Districts (Including All Electoral Districts)",
          label = "tab:ed_boreholes_all",
          covariate.labels = covlabs,
          dep.var.labels = "Number of New Boreholes",
          digits = 2,
          add.lines = list(c("Admin. District Fixed Effects", rep("\\checkmark", 4))),
          omit.stat = c("rsq", "ser", "theta", "aic", "ll"),
          notes.align = "l",
          notes.label = "",
          no.space = TRUE,
          out = "tables/ed_boreholes_all.tex")


# Including only rural electoral districts -------------

ed_rural <- ed %>% filter(urbanprop <= 0.3)

f1 <- formula(cqborehole ~ seg + elf_ed + I(log(popden)) + I(log(urbanprop + 0.001)) +
                log(area_sqkm) + nqborehole98_p10k + admdistrict)

f2 <- formula(cqborehole ~ seg + elf_ed + I(log(popden)) + I(log(urbanprop + 0.001)) +
                log(area_sqkm) + nqborehole98_p10k + avgwinmargpct + coethnic_share +
                prescoethavg + I(log(dist_closest + 1)) + ngowatern_p10k +
                ngoalln_p10k + admdistrict)

# Use function ed_mod to run the formula above for continuous
# segregation and dummied segregation
mods1 <- lapply(c("avediskern", "avediscat3"), function(seg) {
  ed_mod(df = ed_rural, formula = f1, seg = seg, elfco = 0.05)
})

mods2 <- lapply(c("avediskern", "avediscat3"), function(seg) {
  ed_mod(df = ed_rural, formula = f2, seg = seg, elfco = 0.05)
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
          type = "text",
          title = "Segregation and Borehole Investments across Electoral Districts (Rural Electoral Districts Only)",
          label = "tab:ed_boreholes_rural",
          covariate.labels = covlabs,
          dep.var.labels = "Number of New Boreholes",
          digits = 2,
          add.lines = list(c("Admin. District Fixed Effects", rep("\\checkmark", 4))),
          omit.stat = c("rsq", "ser", "theta", "aic", "ll"),
          notes.align = "l",
          notes.label = "",
          no.space = TRUE,
          out = "Replication/tables/ed_boreholes_rural.txt")


# Excluding electoral districts in Machinga and Mangochi -------------

ed_rm <- ed %>% filter(!admdistrict %in% c("Machinga", "Mangochi"))

f1 <- formula(cqborehole ~ seg + elf_ed + I(log(popden)) + I(log(urbanprop + 0.001)) +
                log(area_sqkm) + nqborehole98_p10k + admdistrict)

f2 <- formula(cqborehole ~ seg + elf_ed + I(log(popden)) + I(log(urbanprop + 0.001)) +
                log(area_sqkm) + nqborehole98_p10k + avgwinmargpct + coethnic_share +
                prescoethavg + I(log(dist_closest + 1)) + ngowatern_p10k +
                ngoalln_p10k + admdistrict)

# Use function ed_mod to run the formula above for continuous
# segregation and dummied segregation
mods1 <- lapply(c("avediskern", "avediscat3"), function(seg) {
  ed_mod(df = ed_rm, formula = f1, seg = seg, elfco = 0.05)
})

mods2 <- lapply(c("avediskern", "avediscat3"), function(seg) {
  ed_mod(df = ed_rm, formula = f2, seg = seg, elfco = 0.05)
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
          title = "Segregation and Borehole Investments across Electoral Districts (Excluding Electoral Districts in Machinga and Mangochi)",
          label = "tab:ed_boreholes_rm",
          covariate.labels = covlabs,
          dep.var.labels = "Number of New Boreholes",
          digits = 2,
          add.lines = list(c("Admin. District Fixed Effects", rep("\\checkmark", 4))),
          omit.stat = c("rsq", "ser", "theta", "aic", "ll"),
          notes.align = "l",
          notes.label = "",
          no.space = TRUE,
          out = "tables/ed_boreholes_rm.tex")


# Schools and clinics -------------

f1 <- formula(cqschool ~ seg + elf_ed + I(log(popden)) + I(log(urbanprop + 0.001)) +
                log(area_sqkm) + nqschool98_p10k + admdistrict)

f2 <- formula(cqschool ~ seg + elf_ed + I(log(popden)) + I(log(urbanprop + 0.001)) +
                log(area_sqkm) + nqschool98_p10k + avgwinmargpct + coethnic_share + prescoethavg +
                I(log(dist_closest + 1)) + ngoedun_p10k + ngoalln_p10k + admdistrict)


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
             "Education Aid Projects per 10,000 residents",
             "All Aid Projects per 10,000 residents"
             )

stargazer(mods1[1], mods2[1], mods1[2], mods2[2],
          omit = c("admdistrict"),
          align = T,
          font.size = "small",
          title = "Segregation and School Investments across Electoral Districts",
          label = "tab:ed_schools",
          covariate.labels = covlabs,
          dep.var.labels = "Number of New Schools",
          digits = 2,
          add.lines = list(c("Admin. District Fixed Effects", rep("\\checkmark", 4))),
          omit.stat = c("rsq", "ser", "theta", "aic", "ll"),
          notes.align = "l",
          notes.label = "",
          no.space = TRUE,
          out = "tables/ed_schools.tex")

f1 <- formula(cqclinic ~ seg + elf_ed + I(log(popden)) + I(log(urbanprop + 0.001)) +
                log(area_sqkm) + nqclinic98_p10k + admdistrict)

f2 <- formula(cqclinic ~ seg + elf_ed + I(log(popden)) + I(log(urbanprop + 0.001)) +
                log(area_sqkm) + nqclinic98_p10k + avgwinmargpct+  coethnic_share +
                prescoethavg + I(log(dist_closest + 1)) + ngohealthn_p10k +
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
             "Health Aid Projects per 10,000 residents",
             "All Aid Projects per 10,000 residents"
             )

stargazer(mods1[1], mods2[1], mods1[2], mods2[2],
          omit = c("admdistrict"),
          align = T,
          font.size = "small",
          title = "Segregation and Clinic Investments across Electoral Districts",
          label = "tab:ed_clinics",
          covariate.labels = covlabs,
          dep.var.labels = "Number of New Clinics",
          digits = 2,
          add.lines = list(c("Admin. District Fixed Effects", rep("\\checkmark", 4))),
          omit.stat = c("rsq", "ser", "theta", "aic", "ll"),
          notes.align = "l",
          notes.label = "",
          no.space = TRUE,
          out = "tables/ed_clinics.tex")
