# =========================================================
# data_2020.R
# 2020 US Election Data
# =========================================================

cat("Loading 2020 Election Data...\n")

# ---------------------------------------------------------
# Full dataset
# ---------------------------------------------------------

election2020_full <- read.csv(
  "2020USElection.csv",
  header = TRUE
)

# ---------------------------------------------------------
# Dataset excluding Covid period
# ---------------------------------------------------------

election2020_nocovid <- read.csv(
  "2020USElectionwithoutcovid.csv",
  header = TRUE
)

# ---------------------------------------------------------
# Candidate probabilities
# Use FULL dataset
# ---------------------------------------------------------

Trump_2020 <- as.numeric(
  election2020_full$probTrump
)

Biden_2020 <- as.numeric(
  election2020_full$probBiden
)

# ---------------------------------------------------------
# OU process input
# Use NO-COVID dataset
# ---------------------------------------------------------

Y_2020 <- as.numeric(
  election2020_nocovid$Total
)

# ---------------------------------------------------------
# Time index
# Full timeline
# ---------------------------------------------------------

Sys.setlocale(
  "LC_TIME",
  "English"
)

time_2020 <- as.POSIXct(
  election2020_full$timestampLON,
  format = "%Y/%m/%d %H:%M"
)

# ---------------------------------------------------------
# Covid gap information
# ---------------------------------------------------------

covid_insert_after <- 11593
covid_gap <- nrow(election2020_full) - nrow(election2020_nocovid)

covid_start_full <- covid_insert_after + 1
covid_end_full <- covid_insert_after + covid_gap


# ---------------------------------------------------------
# Mapping:
# No-Covid index -> Full dataset index
# ---------------------------------------------------------

convert_2020_index_to_full <- function(idx) {
  ifelse(
    idx <= covid_insert_after,
    idx,
    idx + covid_gap
  )
}

# ---------------------------------------------------------
# Mapping:
# Full dataset index -> No-Covid index
# ---------------------------------------------------------

convert_2020_full_to_nocovid <- function(idx) {
  
  out <- idx
  
  out[idx > covid_end_full] <-
    idx[idx > covid_end_full] - covid_gap
  
  out
}

# ---------------------------------------------------------
# Basic information
# ---------------------------------------------------------

n_2020_full <- length(time_2020)

n_2020_nocovid <- length(Y_2020)

cat("2020 Full observations:", n_2020_full, "\n")
cat("2020 No-Covid observations:", n_2020_nocovid, "\n")
cat("Covid gap length:", covid_gap, "\n")