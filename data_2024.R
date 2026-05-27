# =========================
# 2024 Election Data
# =========================

b <- read.csv(
  "2024USElection.csv",
  header = TRUE
)

Trump_2024 <- as.numeric(
  b$Donald.Trump.Prob
)

Harris_2024 <- as.numeric(
  b$Kamala.Harris.Prob
)

Y_2024 <- b$Prob

Sys.setlocale("LC_TIME", "English")
time_2024 <- as.POSIXct(
  b$timestampLON,
  format = "%Y/%m/%d %H:%M"
)