# functions.R: A set of custom functions used for the analyses.

## I had to add this in order to make the themes function work below

library(ggplot2)

## This code chunk creates the regression models that we need to use to look for causality

# Function that takes a df, formula, seg variable, and elf cutoff
# and returns a glm (quasi-poisson) object
ed_mod <- function(df, formula, seg, elfco) {
  df <- df %>% plyr::rename(replace = setNames("seg", seg))
  df <- filter(df, elf_ed >= elfco)
  mod <- glm(formula, data = df, family = "quasipoisson")
  return(mod)
}

## This code chunk takes the model created by the function above and creates a matrix with the standard errors

# Function that takes a GLM model and returns (clustered) standard errors or
# the full variance-covariance matrix
ed_se <- function(mod, clustvar = "admdistrict", output = "clustse") {
  fc <- mod$model[, clustvar]
  m <- length(unique(fc))
  k <- length(coef(mod))

  u <- estfun(mod)
  u.clust <- matrix(NA, nrow=m, ncol=k)

  for(j in 1:k){
    u.clust[, j] <- tapply(u[, j], fc, sum)
  }

  vcovm <- vcov(mod)
  cl.vcov <- vcovm %*% ((m / (m-1)) * t(u.clust) %*% (u.clust)) %*% + vcovm

  if(output == "regse") {
    result <- sqrt(diag(vcovm))
  } else if(output == "clustse") {
    result <- sqrt(diag(cl.vcov))
  } else if(output == "regvcov") {
    result <- vcovm
  } else if(output == "clustvcov") {
    result <- cl.vcov
  }

  return(result)
}

## create a function to look at the differences of first differences and predicted values.

# Simulating first differences or predicted values
ed_sim <- function(m, model, vcov, xset_low, xset_high, type = "first-dif") {
  betas <- MASS::mvrnorm(n = m, mu = coef(model), Sigma = vcov)
  pred <- rep(NA, m)

  if (type == "first-dif") {
    for(i in 1:m) {
      pred_low <- exp(xset_low %*% betas[i, ])
      pred_high <- exp(xset_high %*% betas[i, ])
      pred[i] <- pred_high - pred_low
    }
  } else if (type == "pred-val") {
    for(i in 1:m) {
      pred[i] <- exp(xset_low %*% betas[i, ])
    }
  }

  sigquant <- c(0.025, 0.975)

  return(c("mean" = mean(pred), quantile(pred, sigquant)))
}

## Create a graph theme, basically specifying color and size

# A graphing theme for ggplot
theme_smpl <- theme(
    axis.title.x=element_text(size=10, colour='black'),
    axis.title.y=element_text(size=10, angle=90, colour='black'),
    legend.text=element_text(size=9),
    legend.title=element_text(size=9),
    axis.text.x=element_text(size=9),
    axis.text.y=element_text(size=9))
