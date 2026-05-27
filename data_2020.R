# =========================
# 2020 Election Data
# =========================

c <- read.csv(
  "2020USElection.csv",
  header = TRUE
)

d <- read.csv(
  "2020USElectionwithoutcovid.csv",
  header = TRUE
)

Trump_2020 <- as.numeric(
  na.omit(c$probTrump)
)

Biden_2020 <- as.numeric(
  na.omit(c$probBiden)
)

Y_2020 <- d$Total

Sys.setlocale("LC_TIME", "English")
time_2020 <- as.POSIXct(
  c$timestampLON,
  format = "%Y/%m/%d %H:%M"
)