#---------
# Prep
#---------

# Clean
rm(list=ls())

# prepare a seed for simulations
set.seed(2022)

# Set number of runs for our Monte Carlo simulations
MCruns <- 10000

#---------
# Libraries
#---------

# If you don't have a package, this will install it...
if (!require("pacman")) (install.packages("pacman"))
pacman::p_load(
  broom,      # for tidy outputs
  cowplot,    # for plots
  devtools,   # in case you need to install the next package...
  dmetar,     # may need to use devtools::install_github("MathiasHarrer/dmetar")
  effectsize, # for effect sizes
  esc,        # for comparing d values
  googlesheets4, # for loading from googlesheets. Should become obsolete in future version.
  haven,      # for read_dta
  HDInterval,
  Hmisc,
  lazyeval, # for own functions
  lme4,       # for linear mixed effects models
  lmerTest,   # will give significance tests of fixed effects
  lubridate,  # for dates
  magrittr, # For pipe operator %<>%
  meta,     # for metagen
  metafor,  # for meta-analytic averages, metaregressions and forest plots
  metasens, # for limit meta-analysis
  msm,      # for truncated simulations
  patchwork,# for multigraphs
  pracma,
  tidyverse,# For many helpful packages. Like a language within R.
  weights
)
# If a package doesn't load, you may need to update your version of R.

# get rid of summarise messages
options(dplyr.summarise.inform = FALSE)

#---------
# Commands
#---------

# Not in
'%ni%' <- Negate('%in%')

#---------
# Special
#---------

# Just method declare to help get rid of NaNs
# https://stackoverflow.com/questions/52490552/r-convert-nan-to-na
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

#---------
# Functions
#---------

# Get number of unique instances
uniqueN <- function(x) {
  return (length(unique(x)))
}

# Get standard error
se <- function(x) sd(x)/sqrt(length(x))

# Get weighted SD
wtd.sd <- function(x, wt, na.rm = F) {
  sqrt(Hmisc::wtd.var(x, wt, na.rm))
}

# Get pooled SD
getPooledSD <- function(n1, sd1, n2, sd2){
  return(
    abs(
      sqrt(
        (((n1-1)*(sd1^2))+((n2-1)*(sd2^2))) / (n1 + n2 - 2)
      )
    )
  )
}

# Get Cohen's D
getCohenD <- function(m1, m2, pooledSD){
  return(
    (m1 - m2) / pooledSD
  )
}

# Get standard error of Cohen's D
# Done this way in metaregression R textbook
getDSE <- function(d, n1, n2){
  return(
    sqrt(((n1 + n2)/ (n1*n2)) + ((d^2) / (2*(n1 + n2))))
  )
}

# Get share of heterogeneity in multilevel meta-regression
# based on Assink, M., & Wibbelink, C. J. (2016).
# Fitting three-level meta-analytic models in R: A step-by-step tutorial.
metaVar <- function(dat, d_se, model){
  # substitute makes the function like tidyverse functions, you just write
  # the name of the column as a variable.
  
  # Get variance data
  data.variance <- as.data.frame(dat)[,deparse(substitute(d_se))]
  n <- length(data.variance)
  
  # Get variance
  list.inverse.variances <- 1 / (data.variance)
  sum.inverse.variances <- sum(list.inverse.variances)
  squared.sum.inverse.variances <- (sum.inverse.variances) ^ 2
  list.inverse.variances.square <- 1 / (data.variance^2)
  sum.inverse.variances.square <-
    sum(list.inverse.variances.square)
  numerator <- (n - 1) * sum.inverse.variances
  denominator <- squared.sum.inverse.variances -
    sum.inverse.variances.square
  estimated.sampling.variance <- numerator / denominator
  
  # Get I2 values
  I2_denominator <- sum(model$sigma2) + estimated.sampling.variance
  
  returnDf <- NULL
  for (i in 0:length(model$sigma2)) {
    selector <- ((length(model$sigma2) - i) + 1)
    I2_numerator <- ifelse(i == 0,
                           estimated.sampling.variance,
                           model$sigma2[selector])
    
    newDfLine <- data.frame(
      level = i+1,
      name = ifelse(i == 0, "sampling error", model[["s.names"]][selector]),
      T2 = ifelse(i == 0, NA, model$sigma2[selector]),
      I2 = (I2_numerator/I2_denominator)*100
    )
    
    returnDf <- rbind(returnDf, newDfLine)
  }
  
  return(returnDf)
}

# Model comparison is based on:
#  - "Data Analysis: A model Comparison Approach to regression, ANOVA,
# and Beyond" by Judd et al. (2017)
#  - https://github.com/StatQuest/logistic_regression_demo/blob/master/logistic_regression_demo.R

# Custom function to compare to logistic models and get PRD and p
logLik.test.mv <- function (model1, model2){
  
  # Get which model is the reduced (C) or full (A)
  if(length(coef(model1)) <= length(coef(model2))) {
    mC <- model1
    mA <- model2
  } else {
    mC <- model2
    mA <- model1
  }
  
  df <- length(coef(mA)) - length(coef(mC))
  LLC <- -2*logLik(mC)
  LLA <- -2*logLik(mA)
  diffLL <- abs(LLC - LLA) # abs to deal with rare case of positive logLik
  PRD <- abs(diffLL / LLC)
  
  critChi <- qchisq(1-.05, df = df)
  p <- pchisq(diffLL, df=df, lower.tail=FALSE)
  
  outputTable <- data.frame(
    PRD=round.c(PRD, digits = 3),
    diffLL=round.c(diffLL, digits = 2),
    df=df,
    critChi=round.c(critChi, digits = 2),
    p=round.c(p, digits = 3)
  )
  
  return(
    list(
      modelC = mC$formula.mods,
      modelA = mA$formula.mods,
      outputTable = outputTable
    )
  )
}

# Function to get SE from mean and CI
getSEfromCI <- function(CIlow, CIup){
  return((CIup - CIlow)/3.92)
}

# Round number to string, fills in missing digits (decimals) with 0s
round.c <- function(x, digits = 2) {
  roundedElement <- round(x, digits)
  roundedElement.char <- as.character(roundedElement)
  
  # If digits are more than 0
  if(digits > 0) {
    # If no decimals
    if(roundedElement %% 1 == 0) {
      roundedElement.char <- paste0(roundedElement.char, ".")
      for (i in 1:digits) {
        roundedElement.char <- paste0(roundedElement.char, "0")
      }
    } else {
      ndecimals <- nchar(strsplit(roundedElement.char, ".", fixed = T)[[1]][2])
      diffToDigits <- digits - ndecimals
      
      if(diffToDigits > 0) {
        for (i in 1:diffToDigits) {
          roundedElement.char <- paste0(roundedElement.char, "0")
        }
      }
    }
  }
  
  return(roundedElement.char)
}

# Function to report a line from a multi-level meta-regression #
reportMVReg <- function(model, index){
  
  # Dealing with special bits of the p value
  pvalue <- round.c(model$pval[index], 3)
  if(as.numeric(pvalue) == 1){
    pvalue <- " = 1"
  } else {
    prefix <- ifelse(as.numeric(pvalue) == 0, " < .", " = .")
    pvalue <- ifelse(as.numeric(pvalue) == 0, "0.001", pvalue)
    pvalue <- paste0(prefix, strsplit(pvalue, ".", fixed = T)[[1]][2])
  }
  
  return (
    paste0(
      "b = ", round.c(coef(model)[index], 2),
      ", SE = ", round.c(model$se[index], 2),
      ", t(", model$ddf[index], ") = ", round.c(model$zval[index], 2),
      ", p", pvalue,
      ", 95% CI[", round.c(model$ci.lb[index], 2), ", ",
      round.c(model$ci.ub[index], 2), "]."
    )
  )
}

# Function to report a distribution #
MMSDCI.vec <- function(x, ci = .95){
  
  # Get the upper and lower quantiles
  CI.low <- (1 - ci)/2
  CI.upp <- ci + CI.low
  
  # Return the summary
  return(
    paste0(
      "M = ", round.c(mean(x, na.rm=T), 2),
      ", SD = ", round.c(sd(x, na.rm=T), 2),
      ", Median = ", round.c(median(x, na.rm=T), 2),
      ", (95% CI: ", round.c(quantile(x, probs = c(CI.low, CI.upp))[[1]], 2), ", ",
      round.c(quantile(x, probs = c(CI.low, CI.upp))[[2]], 2), ")"
    )
  )
}

MMSDCI <- function(x, ci = .95){
  
  # If it is vector, just report one line
  if (is.vector(x)){
    report <- MMSDCI.vec(x, ci)
  } else {
    # Filter out non-numeric columns
    x <- x %>% select_if(is.numeric)
    
    # Otherwise, report a line for every column
    report <- lapply(x, MMSDCI.vec)
  }
  
  return(report)
}

# Make into a percentage string
makeIntoPercString <- function(x, digits = 2) {
  return(paste0(round.c(x*100, digits), "%"))
}

# Combine PE and SIM
combinePEandSIM <- function(pe, sim, ci = .95){
  
  # Get the upper and lower quantiles
  CI.low <- (1 - ci)/2
  CI.upp <- ci + CI.low
  
  # Prepare the new dataframe where each row is a variable
  # And we will have a pe, ci, and combined, column
  df <- pe %>% pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "pe"
  )
  
  # Get the CI for each variable
  df$ci <- NA
  for (variable in df$variable) {
    
    newCI <- paste0(
      "(95% CI: ", 
      round.c(quantile(as.data.frame(sim)[,variable], 
                       probs = c(CI.low, CI.upp))[[1]], 2), ", ",
      round.c(quantile(as.data.frame(sim)[,variable], 
                       probs = c(CI.low, CI.upp))[[2]], 2), ")"
    )
    
    df$ci[which(df$variable == variable)] <- newCI
  }
  
  # Copy into one string, for ease fo copy pasting
  df <- df %>% rowwise() %>% 
    mutate(combined = paste0(round.c(pe, 2), " ", ci)) %>% 
    ungroup()
  
  return(df)
}
