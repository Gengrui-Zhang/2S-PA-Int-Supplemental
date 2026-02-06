# =============================================================================
# 2S-PA-Int Demonstration with Synthetic Data of PIRLS 2021
#
# Note (to include in methods / script comments):
# For the main empirical results in the manuscript, we used the Croatia
# sample from PIRLS 2021 (Mullis et al., 2023). Because PIRLS 2021 public
# data are externally hosted and distributed, and to comply with journal TOP
# Level 2 transparency requirements while preserving data privacy and
# licensing conditions, this script generates a synthetic dataset using
# the synthpop package and demonstrates the key analysis using that synthetic
# data. Users can obtain the original PIRLS public-use data from the
# official PIRLS 2021 site if desired.
# =============================================================================

# -----------------------------------------------------------------------------
# Load required libraries
# -----------------------------------------------------------------------------
library(here)
library(haven)
library(dplyr)
library(synthpop)
library(lavaan)
library(psych)

# -----------------------------------------------------------------------------
# 1) Load the original PIRLS SPSS dataset (public-use variables)
# -----------------------------------------------------------------------------
# Read the original public data
PIRLS2021 <- read_sav(here("Paper", "PIRLS2021", "ASGHRVA5.sav"))

# -----------------------------------------------------------------------------
# 2) Select relevant variables and drop missing
# -----------------------------------------------------------------------------
PIRLS_Data <- PIRLS2021 %>%
  select(ASBR07A:ASBR07H, ASRREA01:ASRREA05) %>%
  drop_na()

# -----------------------------------------------------------------------------
# 3) Remove SPSS value labels but preserve numeric codes
# -----------------------------------------------------------------------------
PIRLS_Data_clean <- PIRLS_Data %>%
  mutate(across(where(haven::is.labelled), haven::zap_labels))

# -----------------------------------------------------------------------------
# 4) Convert small-level numeric variables to ordered factors
# -----------------------------------------------------------------------------
PIRLS_Data_clean <- PIRLS_Data_clean %>%
  mutate(
    across(
      c(ASBR07A, ASBR07B, ASBR07C, ASBR07D,
        ASBR07E, ASBR07F, ASBR07G, ASBR07H),
      ~ factor(., ordered = TRUE)
    )
  )

# -----------------------------------------------------------------------------
# 5) Generate one synthetic dataset
# -----------------------------------------------------------------------------
# Demonstration-only: Use synthpop to generate one synthetic dataset
syn_pirls <- syn(
  PIRLS_Data_clean,
  m = 1
)

# -----------------------------------------------------------------------------
# 6) Analysis demonstration: 2S-PA-Int Model on synthetic data
# -----------------------------------------------------------------------------
# 6a) Recode variables to align with original analysis coding
# Read the synthetic data
PIRLS_synthetic <- read_csv(here("Paper", "PIRLS2021", "PIRLS_synthetic.csv"))

syn_data <- PIRLS_synthetic %>%
  mutate(across(
    c(ASBR07A, ASBR07B, ASBR07D, ASBR07E, ASBR07F, ASBR07G, ASBR07H),
    ~ recode(.x, `1` = 4, `2` = 3, `3` = 2, `4` = 1)
  ))

# 6b) Prepare measurement constructs
IM <- syn_data %>%
  select(ASBR07C, ASBR07D, ASBR07E) %>%
  rename(IM1 = ASBR07C, IM2 = ASBR07D, IM3 = ASBR07E)

EM <- syn_data %>%
  select(ASBR07A, ASBR07B, ASBR07F) %>%
  rename(EM1 = ASBR07A, EM2 = ASBR07B, EM3 = ASBR07F)

DV <- syn_data$ASRREA01
final_df <- cbind(IM, EM, DV)

# 6c) Obtain factor scores
mod_cfa <- '
  IM =~ IM1 + IM2 + IM3
  EM =~ EM1 + EM2 + EM3
'

fs_dat <- get_fs(final_df, model = mod_cfa,
                 method = "Bartlett", std.lv = TRUE)

# Create interaction latent product
fs_dat$fs_int <- fs_dat$fs_IM * fs_dat$fs_EM
fs_dat$fs_int <- fs_dat$fs_int - mean(fs_dat$fs_int)
fs_dat$fs_int_se <- sqrt(
  fs_dat$fs_IM_se^2 +
    fs_dat$fs_EM_se^2 +
    fs_dat$fs_IM_se^2 * fs_dat$fs_EM_se^2
)

fs_dat$DV <- DV

# 6d) Specify 2S-PA-Int SEM model with error constraints
sem_model <- "
  IM =~ 1*fs_IM
  EM =~ 1*fs_EM
  lInt =~ 1*fs_int

  fs_IM ~~ ev_IM*fs_IM
  fs_EM ~~ ev_EM*fs_EM
  fs_int ~~ ev_int*fs_int

  DV ~ b1*IM + b2*EM + b3*lInt

  IM ~~ v1*IM
  EM ~~ v2*EM
  lInt ~~ v3*lInt

  IM ~~ cov_12*EM
  IM ~~ cov_13*lInt
  EM ~~ cov_23*lInt

  DV ~~ dist_y*DV

  var_y := dist_y +
            (b1^2 * v1 + b2^2 * v2 + b3^2 * v3 +
             2*b1*b2*cov_12 + 2*b1*b3*cov_13 + 2*b2*b3*cov_23)

  beta1 := b1*sqrt(v1)/sqrt(var_y)
  beta2 := b2*sqrt(v2)/sqrt(var_y)
  beta3 := b3*sqrt(v1)*sqrt(v2)/sqrt(var_y)
"

sem_model <- gsub("ev_IM", fs_dat$fs_IM_se[1]^2, sem_model)
sem_model <- gsub("ev_EM", fs_dat$fs_EM_se[1]^2, sem_model)
sem_model <- gsub("ev_int", fs_dat$fs_int_se[1]^2, sem_model)

# Fit the model
fit_tspa <- sem(sem_model, data = fs_dat)

# 7e) Output standardized parameter estimates
est <- parameterEstimates(fit_tspa, standardized = TRUE)

results <- est[est$label %in% c("beta1", "beta2", "beta3"), ]
print(results)

cat("\n2S-PA-Int analysis on synthetic data complete.\n")

# =============================================================================
# End of script
# =============================================================================
