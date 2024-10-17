# Illustrative Example of 2S-PA-Int
# Load packages
library(haven)
library(here)
library(lavaan)
library(psych)

# Read Data
PIRLS2021 <- read_sav(here("PIRLS2021", "ASGHRVA5.sav"))

# Select relevant variables
PIRLS_Data <- PIRLS2021 %>%
  select(ASBR07A:ASBR07H, ASRREA01:ASRREA05) %>%
  drop_na()

# Recode variables
PIRLS_Data <- PIRLS_Data %>%
  mutate(across(
    .cols = c(ASBR07A:ASBR07B, ASBR07D:ASBR07H),
    .fns = ~ recode(as.numeric(.x), `1` = 4, `2` = 3, `3` = 2, `4` = 1),
    .names = "{.col}_recode"
  )) %>%
  drop_na()

# Measurement Model
IM <- PIRLS_Data %>% 
  select(ASBR07C, ASBR07D_recode, ASBR07E_recode) %>%
  rename(
    IM1 = ASBR07C,
    IM2 = ASBR07D_recode,
    IM3 = ASBR07E_recode
  )

EM <- PIRLS_Data %>% 
  select(ASBR07A_recode, ASBR07B_recode, ASBR07F_recode) %>%
  rename(
    EM1 = ASBR07A_recode,
    EM2 = ASBR07B_recode,
    EM3 = ASBR07F_recode
  )

model <- "IM =~ IM1 + IM2 + IM3
          EM =~ EM1 + EM2 + EM3"

# Add reading performance score
DV <- "ASRREA01"
final_df <- cbind(IM, EM, PIRLS_Data[DV])

# Assess measurement model fit
fit <- fitmeasures(sem(model, final_df))

# Fit 2S-PA-Int Model
# Stage 1: Obtain factor score estimates with standard errors of measurement
mod_cfa <- '
          IM =~ IM1 + IM2 + IM3
          EM =~ EM1 + EM2 + EM3
'
# Obtain factor scores
fs_dat <- get_fs(final_df, model = mod_cfa, method = "Bartlett", std.lv = TRUE)
# Obtain factor product
fs_dat$fs_int <- fs_dat$fs_IM * fs_dat$fs_EM
# Centering the product variable
fs_dat$fs_int <- fs_dat$fs_int - mean(fs_dat$fs_int)
fs_dat$fs_int_se <- sqrt(1 * fs_dat$fs_IM_se^2 + 1 * fs_dat$fs_EM_se^2 + fs_dat$fs_IM_se^2*fs_dat$fs_EM_se^2)

# Stage 2: Fit structural model with interaction
mod_tspa <- "
            # latent variables (indicated by factor scores)
            IM =~ 1 * fs_IM
            EM =~ 1 * fs_EM
            lInt =~ 1 * fs_int

            # constrain the errors
            fs_IM ~~ ev_IM * fs_IM
            fs_EM ~~ ev_EM * fs_EM
            fs_int ~~ ev_int * fs_int
            
            # regressions
            DV ~ b1*IM + b2*EM + b3*lInt

            # Variance
            IM ~~ v1*IM
            EM ~~ v2*EM
            lInt ~~ v3*lInt
            
            # Covariance
            IM ~~ cov_12*EM
            IM ~~ cov_13*lInt
            EM ~~ cov_23*lInt
          
            # Disturbance
    			  DV ~~ dist_y*DV
            var_y :=  dist_y + (b1^2 * v1 + b2^2 * v2 + b3^2 * v3 + 
                              2 * b1 * b2 * cov_12 + 
                              2 * b1 * b3 * cov_13 + 
                              2 * b2 * b3 * cov_23)
              
            # Define Standardized Coefficients
            beta1 := b1*sqrt(v1)/sqrt(var_y)
            beta2 := b2*sqrt(v2)/sqrt(var_y)
            beta3 := b3*sqrt(v1)*sqrt(v2)/sqrt(var_y)"

# Replace with error-variance constraints
mod_tspa <- gsub("ev_IM", replacement = fs_dat$fs_IM_se[1]^2, x = mod_tspa)
mod_tspa <- gsub("ev_EM", replacement = fs_dat$fs_EM_se[1]^2, x = mod_tspa)
mod_tspa <- gsub("ev_int", replacement = fs_dat$fs_int_se[1]^2, x = mod_tspa)

# Append reading performance score
fs_dat$DV <- PIRLS_Data[DV]
# Fit 2S-PA-Int Model
fit_tspa <- sem(mod_tspa,
                data = fs_dat)

# Extract parameters
beta1_tspa <- c(beta1_tspa, coef(fit_tspa, type = "user")["beta1"])
beta2_tspa <- c(beta2_tspa, coef(fit_tspa, type = "user")["beta2"])
beta3_tspa <- c(beta3_tspa, coef(fit_tspa, type = "user")["beta3"])
se1_tspa <- c(se1_tspa, sqrt(vcov(fit_tspa, type = "user")["beta1", "beta1"]))
se2_tspa <- c(se2_tspa, sqrt(vcov(fit_tspa, type = "user")["beta2", "beta2"]))
se3_tspa <- c(se3_tspa, sqrt(vcov(fit_tspa, type = "user")["beta3", "beta3"]))
