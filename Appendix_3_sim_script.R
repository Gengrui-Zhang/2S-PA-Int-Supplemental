# Load pacakges
library(semTools)
library(lavaan)
library(SimDesign)
library(mnormt)
library(dplyr)
library(tidyverse)
library(semlrtp)
library(here)

# Source Function
r_scripts <- list.files(here("R"), pattern = "\\.R$", full.names = TRUE)
lapply(r_scripts, source)

# ========================================= Simulation Conditions ========================================= #

DESIGNFACTOR <- createDesign(
  N = c(100, 250, 500),
  cor_xm = c(0, 0.3, 0.6), # correlation between latent x and m / error variance of Y
  rel = c(0.7, 0.8, 0.9),
  gamma_xm = c(0, 0.3) # Two levels of the interaction effect
)

FIXED_PARAMETER <- list(model = '
                                  # Measurement Model
                                    X =~ x1 + x2 + x3
                                    M =~ m1 + m2 + m3
                                  # Structural Model
                                    Y ~ b1*X + b2*M + b3*X:M
                                  # Define Standardized Coefficients
                                    X ~~ v1*X
                                    M ~~ v2*M
                                    beta1 := b1*sqrt(v1)
                                    beta2 := b2*sqrt(v2)
                                    beta3 := b3*sqrt(v1)*sqrt(v2)
                                  ',
                        beta1 = 1,  # fixed
                        beta2 = 0.9,  # fixed
                        beta3 = 0.75,  # three conditions
                        mu_x = 0, # latent mean of x: fixed at 0
                        mu_m = 0, # latent mean of m: fixed at 0,
                        gamma_x = 0.3,
                        gamma_m = 0.3
)

# ========================================= Data Generation ========================================= #
# Helper Function
generate_sem_data <- function(N, model, Alpha, Phi, Lambda, Gamma, Theta, SD_y) {
  # Generate scores for observed items: x1 - x3, m1 - m3
  eta_scores <- rmnorm(N, mean = Alpha, varcov = Phi) # Simulate latent scores
  # Generate eta scores with interaction
  eta_scores_int <- cbind(eta_scores, eta_scores[,1]*eta_scores[,2])
  delta <- rmnorm(N, mean = rep(0, length(diag(Theta))), varcov = Theta) # Errors/Residuals of Indicators
  item_scores <- tcrossprod(eta_scores, Lambda) + delta # Item/Indicator scores

  y_scores <- tcrossprod(eta_scores_int, Gamma) + rnorm(N, 0, SD_y) # Y scores / DV scores

  # Parsing formula
  indicator_terms <- unlist(strsplit(gsub(" ", "",
                                          unlist(strsplit(model, "\n"))[grep("=~",
                                                                             unlist(strsplit(model, "\n")))]),
                                     split = "=~"))
  indicator_vars <- indicator_terms[grepl("+", indicator_terms, fixed = TRUE) == "FALSE"]
  indicator_num <- as.vector(unlist(lapply(strsplit(indicator_terms[grepl("+", indicator_terms, fixed = TRUE) == "TRUE"],
                                                    split = "+", fixed = TRUE), length)))

  df <- as.data.frame(cbind(item_scores, y_scores))
  names <- list()
  for (n in seq(length(indicator_vars))) {
    names[[n]] <- tolower(paste0(indicator_vars[n], seq(indicator_num[n])))
  }
  colnames(df) <- c(unlist(names), "Y")
  return(df)
}

generate_dat <- function (condition, fixed_objects = NULL) {
  N <- condition$N # Sample size
  mu_x <- fixed_objects$mu_x
  mu_m <- fixed_objects$mu_m
  beta1 <- fixed_objects$beta1 # beta 1: fixed at 1
  beta2 <- fixed_objects$beta2 # beta 2: fixed at 0.9
  beta3 <- fixed_objects$beta3 # beta 3: fixed at 0.75
  cor_xm <- condition$cor_xm # latent correlation: varied
  gamma_x <- fixed_objects$gamma_x
  gamma_m <- fixed_objects$gamma_m
  gamma_xm <- condition$gamma_xm
  rel = condition$rel

  # Compute disturbance variance
  sd_y <- sqrt(1 - (gamma_x^2 + gamma_m^2 + gamma_xm^2 + 2*gamma_x*gamma_m*cor_xm))

  # Compute error variance
  sum_error <- sum(c(beta1, beta2, beta3))^2*(1 - rel)/rel
  err_var <- sum_error*c(0.44, 0.33, 0.23)

  # Simulate SEM data
  Alpha <- c(mu_x, mu_m) # Latent means
  Phi <- matrix(c(1, cor_xm,
                  cor_xm, 1), nrow = 2) # latent var/cov
  Lambda <- cbind(c(beta1, beta2, beta3, rep(0, 3)),
                  c(rep(0, 3), beta1, beta2, beta3)) # factor loadings
  Theta <- diag(err_var,
                nrow = 6)
  Gamma <- matrix(c(gamma_x, gamma_m, gamma_xm), nrow = 1)
  SD_y <- sd_y

  generate_sem_data(N,
                    model = fixed_objects$model,
                    Alpha = Alpha,
                    Phi = Phi,
                    Lambda = Lambda,
                    Theta = Theta,
                    Gamma = Gamma,
                    SD_y = SD_y
  )
}

# ========================================= Data Analysis ========================================= #
# Moderated Multiple Regression
analyze_mmr <- function(condition, dat, fixed_objects = NULL) {
  
  # Initialize a local warning counter for this replication
  local_warning_counter <- 0
  
  # Fit the model and capture warnings
  result <- withCallingHandlers(
    {
      # Fit using moderated multiple regression
      dat_mmr <- dat %>%
        mutate(x_c = rowSums(dat[c("x1", "x2", "x3")]) - mean(rowSums(dat[c("x1", "x2", "x3")])),
               m_c = rowSums(dat[c("m1", "m2", "m3")]) - mean(rowSums(dat[c("m1", "m2", "m3")])),
               xm = x_c * m_c)
      
      # Fit the model
      fit_mmr <- sem(model = "Y ~ c0*x_c + c1*m_c + c2*xm", 
                     data = dat_mmr,
                     bounds = TRUE)
      
      # Extract parameters
      std_col <- standardizedSolution(fit_mmr)
      est <- std_col[std_col$label == "c2", "est.std"]
      se <- std_col[std_col$label == "c2", "se"]
      usd_col <- parameterEstimates(fit_mmr, standardized = FALSE)
      est_usd <- usd_col[usd_col$label == "c2", "est"]
      se_usd <- usd_col[usd_col$label == "c2", "se"]
      lrtp_col <- lrtp(fit_mmr)
      lrtp_ci <- as.data.frame(lrtp_col) %>%
        filter(op == "~" & label == "c2") %>%
        select(ci.lower, ci.upper)
      
      # Create the output vector
      out <- c(est, se, est_usd, se_usd, lrtp_ci$ci.lower, lrtp_ci$ci.upper, local_warning_counter)  
      names(out) <- c("est", "se_std", "est_usd", "se_usd", "lrtp_lower", "lrtp_upper", "warnings_count")
      
      # Check if any NA exists in the output
      if (anyNA(out)) {
        stop("The model did not obtain SE")
      }
      
      return(out)
    },
    warning = function(w) {
      # Increment the local warning counter
      local_warning_counter <<- local_warning_counter + 1
    }
  )
  
  return(result)
}
# Matched-Pair Unconstrained Product Indicator
analyze_upi <- function(condition, dat, fixed_objects = NULL) {
  
  # Initialize a local warning counter for this replication
  local_warning_counter <- 0
  
  # Fit the model and capture warnings
  result <- withCallingHandlers(
    {
      fit_upi <- upi(model = fixed_objects$model, 
                     data = dat, 
                     mode = "match",
                     bounds = TRUE) 
      
      # Extract parameters
      est <- coef(fit_upi, type = "user")["beta3"]
      se <- sqrt(vcov(fit_upi, type = "user")["beta3", "beta3"])
      est_usd <- coef(fit_upi)["b3"]
      se_usd <- sqrt(vcov(fit_upi)["b3", "b3"])
      lrtp_col <- lrtp(fit_upi)
      lrtp_ci <- as.data.frame(lrtp_col) %>%
        filter(op == "~" & label == "b3") %>%
        select(ci.lower, ci.upper)
      
      # Create the output vector
      out <- c(est, se, est_usd, se_usd, lrtp_ci$ci.lower, lrtp_ci$ci.upper, local_warning_counter)  
      names(out) <- c("est", "se_std", "est_usd", "se_usd", "lrtp_lower", "lrtp_upper", "warnings_count")
      
      # Check for NA values in the output
      if (anyNA(out)) {
        stop("The model did not obtain SE")
      }
      
      return(out)
    },
    warning = function(w) {
      # Increment the local warning counter
      local_warning_counter <<- local_warning_counter + 1
    }
  )
  
  return(result)
}
# Reliability-Adjusted Product Indicator
analyze_rapi <- function(condition, dat, fixed_objects = NULL) {
  
  # Initialize a local warning counter for this replication
  local_warning_counter <- 0
  
  # Fit the model and capture warnings
  result <- withCallingHandlers(
    {
      # Fit the model
      fit_rapi <- rapi(model = fixed_objects$model,
                       data = dat,
                       bounds = TRUE)
      
      # Extract parameters
      est <- coef(fit_rapi, type = "user")["beta3"]
      se <- sqrt(vcov(fit_rapi, type = "user")["beta3", "beta3"])
      est_usd <- coef(fit_rapi)["b3"]
      se_usd <- sqrt(vcov(fit_rapi)["b3", "b3"])
      lrtp_col <- lrtp(fit_rapi)
      lrtp_ci <- as.data.frame(lrtp_col) %>%
        filter(op == "~" & label == "b3") %>%
        select(ci.lower, ci.upper)
      
      # Create the output vector
      out <- c(est, se, est_usd, se_usd, lrtp_ci$ci.lower, lrtp_ci$ci.upper, local_warning_counter)
      names(out) <- c("est", "se_std", "est_usd", "se_usd", "lrtp_lower", "lrtp_upper", "warnings_count")
      
      # Check for NA values in the output
      if (anyNA(out)) {
        stop("The model did not obtain SE")
      }
      
      return(out)
    },
    warning = function(w) {
      # Increment the local warning counter
      local_warning_counter <<- local_warning_counter + 1
    }
  )
  
  return(result)
}
# Two-Stage Path Analysis with Interaction
analyze_tspa <- function(condition, dat, fixed_objects = NULL) {
  
  # Initialize a local warning counter for this replication
  local_warning_counter <- 0
  
  # Fit the model and capture warnings
  result <- withCallingHandlers(
    {
      # Fit using two-stage path analysis
      fs_dat <- get_fs(dat,
                       model = '
                         X =~ x1 + x2 + x3
                         M =~ m1 + m2 + m3
                       ',
                       method = "Bartlett",
                       std.lv = TRUE,
                       bounds = TRUE)
      
      Y <- dat$Y
      fs_dat <- cbind(fs_dat, Y)
      
      fit_tspa <- tspa(model = "Y ~ b1*X + b2*M + b3*X:M
                                beta1 := b1 * sqrt(v1)
                                beta2 := b2 * sqrt(v2)
                                beta3 := b3 * sqrt(v1) * sqrt(v2)",
                       data = fs_dat,
                       se = list(X = fs_dat$fs_X_se[1], M = fs_dat$fs_M_se[1]),
                       bounds = TRUE) 
      
      # Extract parameters
      est <- coef(fit_tspa, type = "user")["beta3"]
      se <- sqrt(vcov(fit_tspa, type = "user")["beta3", "beta3"])
      est_usd <- coef(fit_tspa)["b3"]
      se_usd <- sqrt(vcov(fit_tspa)["b3", "b3"])
      lrtp_col <- lrtp(fit_tspa)
      lrtp_ci <- as.data.frame(lrtp_col) %>%
        filter(op == "~" & label == "b3") %>%
        select(ci.lower, ci.upper)
      
      # Create the output vector
      out <- c(est, se, est_usd, se_usd, lrtp_ci$ci.lower, lrtp_ci$ci.upper, local_warning_counter)
      names(out) <- c("est", "se_std", "est_usd", "se_usd", "lrtp_lower", "lrtp_upper", "warnings_count")
      
      # Check for NA values in the output
      if (anyNA(out)) {
        stop("The model did not obtain SE")
      }
      
      return(out)
    },
    warning = function(w) {
      # Increment the local warning counter
      local_warning_counter <<- local_warning_counter + 1
    }
  )
  
  return(result)
}

# ========================================= Results Summary ========================================= #
# Helper function: robust bias
robust_bias <- function(est, se, par, trim = 0, type = NULL) {
  output <- numeric(ncol(est))
  for (i in seq_len(ncol(est))) {
    if (type == "raw") {
      output[i] <- mean((est[,i] - par), na.rm = TRUE)
    } else if (type == "standardized") {
      output[i] <- (mean(est[,i], na.rm = TRUE) - par)/sd(est[,i], na.rm = TRUE)
    } else if (type == "trim") {
      output[i] <- mean(est[,i], trim = trim, na.rm = TRUE) - par
    } else if (type == "median") {
      output[i] <- (median(est[,i], na.rm = TRUE) - par) / mad(est[,i], na.rm = TRUE)
    } else {
      output[i] <- (mean(est[,i], trim = trim, na.rm = TRUE) - par) / sd(est[,i], na.rm = TRUE)
    }
  }
  names(output) <- colnames(est)
  return(output)
}

# Helper function: relative SE bias
rse_bias <- function(est, se, trim = 0, type = "raw") {
  if (type == "raw") {
    se_mean <- apply(se, 2, mean, na.rm = T)
    se_sd <- apply(est, 2L, sd, na.rm = T)
    rse_bias <- se_mean / se_sd - 1
  } else if (type == "median") {
    se_median <- apply(se, 2, median, na.rm = TRUE)
    mad_sd <- apply(est, 2, function(x) mad(x, na.rm = TRUE))
    rse_bias <- se_median / mad_sd - 1
  } else if (type == "trim") {
    se_mean <- apply(se, 2, mean, trim = trim, na.rm = TRUE)
    se_sd <- apply(est, 2L, sd, na.rm = T)
    rse_bias <- se_mean / se_sd - 1
  }
  return(rse_bias)
}

# Helper function: detecting outliers for SE
outlier_se <- function(se) {
  results <- c()
  for(column in colnames(se)) {
    # Calculate Q1, Q3, and IQR
    Q1 <- quantile(se[,column], 0.25, na.rm = TRUE)
    Q3 <- quantile(se[,column], 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    # Determine outliers
    lower_bound <- (Q1 - 1.5 * IQR)
    upper_bound <- (Q3 + 1.5 * IQR)
    outliers <- se[,column][se[,column] < lower_bound | se[,column] > upper_bound]
    # Calculate the percentage of outliers
    percentage <- length(outliers) / sum(!is.na(se[,column])) * 100
    results[column] <- percentage
  }
  return(results)
}

# Helper function for calculating coverage rate, Type I error rate, and power
ci_stats <- function(est, se, par, stats_type, lrt_lo = NULL, lrt_up = NULL) {

  # Calculate the confidence intervals (usd)
  lo_95 <- est - qnorm(.975) * se
  up_95 <- est + qnorm(.975) * se
  ci_est <- vector("list", length = ncol(est))
  names(ci_est) <- colnames(est)
  
  # Construct confidence intervals for each method
  for (i in seq_len(ncol(est))) {
    ci_est[[i]] <- cbind(lo_95[,i], up_95[,i])
  }
  
  # Extract LRT CIs
  if (!is.null(lrt_lo) && !is.null(lrt_up)) {
    ci_lrt <- vector("list", length = ncol(est))
    names(ci_lrt) <- colnames(est)
    
    for (i in seq_len(ncol(est))) {
      ci_lrt[[i]] <- cbind(lrt_lo[,i], lrt_up[,i])
    }
  }
  
  # Determine which statistic to calculate
  if (stats_type == "Coverage") {
    return(sapply(ci_est, function(ci) mean(ci[,1] <= par & ci[,2] >= par)))
  } else if (stats_type == "TypeI") {
    return(sapply(ci_est, function(ci) mean(ci[,1] > 0 | ci[,2] < 0)))
  } else if (stats_type == "Lrt_TypeI") {
    return(sapply(ci_lrt, function(ci) mean(ci[,1] > 0 | ci[,2] < 0)))
  } else if (stats_type == "Power") {
    return(sapply(ci_est, function(ci) (1 - mean(ci[,1] < 0 & ci[,2] > 0))))
  } else if (stats_type == "Lrt_Power") {
    return(sapply(ci_lrt, function(ci) (1 - mean(ci[,1] < 0 & ci[,2] > 0))))
  } else {
    stop("Invalid stats_type specified. Please choose from 'Coverage', 'TypeI', or 'Power'.")
  }
}

# Helper function for warning sum
warning_sum <- function(count) {
  apply(count, 2, sum, na.rm = TRUE)
}

# Evaluation Function
evaluate_res <- function (condition, results, fixed_objects = NULL) {

  # Population parameter
  pop_par <- condition$gamma_xm

  # Parameter estimates
  est_std <- results[, grep(".est$", colnames(results))]
  est_usd <- results[, grep(".est_usd$", colnames(results))]
  se_std <- results[, grep(".se_std$", colnames(results))]
  se_usd <- results[, grep(".se_usd$", colnames(results))]
  lrtp_lower <- results[, grep(".lrtp_lower$", colnames(results))]
  lrtp_upper <- results[, grep(".lrtp_upper$", colnames(results))]
  warnings <- results[, grep(".warnings_count$", colnames(results))]

  c(raw_bias = robust_bias(est_std,
                           se_std,
                           pop_par,
                           type = "raw"),
    std_bias = robust_bias(est_std,
                           se_std,
                           pop_par,
                           type = "standardized"),
    trim_bias = robust_bias(est_std,
                            se_std,
                            pop_par,
                            trim = 0.2,
                            type = "trim"), # 20% trimmed mean
    stdMed_bias = robust_bias(est_std,
                              se_std,
                              pop_par,
                              type = "median"),
    raw_rse_bias = rse_bias(est_std,
                            se_std,
                            type = "raw"),
    stdMed_rse_bias = rse_bias(est_std,
                               se_std,
                               type = "median"),
    trim_rse_bias = rse_bias(est_std,
                             se_std,
                             trim = 0.2,
                             type = "trim"),
    outlier_se = outlier_se(se_std),
    coverage_usd = ci_stats(est_usd, 
                        se_usd,
                        pop_par, 
                        "Coverage"),
    coverage_std = ci_stats(est_std, 
                            se_std,
                            pop_par, 
                            "Coverage"),
    type1_usd = ci_stats(est_usd, 
                         se_usd,
                         pop_par, 
                         "TypeI"),
    type1_std = ci_stats(est_std, 
                         se_std,
                         pop_par, 
                         "TypeI"),
    type1_lrt = ci_stats(est_usd, 
                         se_usd,
                         pop_par, 
                         "Lrt_TypeI",
                         lrt_lo = lrtp_lower,
                         lrt_up = lrtp_upper),
    power_usd = ci_stats(est_usd, 
                         se_usd,
                         pop_par, 
                         "Power"),
    power_std = ci_stats(est_std, 
                         se_std,
                         pop_par, 
                         "Power"),
    power_lrt = ci_stats(est_usd, 
                         se_usd,
                         pop_par, 
                         "Lrt_Power",
                         lrt_lo = lrtp_lower,
                         lrt_up = lrtp_upper),
    rmse = RMSE(na.omit(est_std),
                parameter = pop_par),
    warning_total = warning_sum(warnings)
  )
}

# ========================================= Run Experiment ========================================= #
res <- runSimulation(design = DESIGNFACTOR,
              replications = 2000,
              generate = generate_dat,
              analyse = list(mmr = analyze_mmr,
                             upi = analyze_upi,
                             rapi = analyze_rapi,
                             tspa = analyze_tspa),
              summarise = evaluate_res,
              fixed_objects = FIXED_PARAMETER,
              seed = rep(61543, nrow(DESIGNFACTOR)),
              packages = "lavaan", 
              filename = "continuous_boundMLR_09252024",
              parallel = TRUE,
              ncores = 30,
              save = TRUE,
              save_results = TRUE)
