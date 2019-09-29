# Code originally from: https://raw.githubusercontent.com/M4Competition/M4-methods/master/Benchmarks%20and%20Evaluation.R
# Modifications and clean up by Gary Mulder mulder.g@unic.ac.cy
# Version 0.1 - 2019-09-24

library(forecast) #Requires v8.2

smape <- function(fcast, out_sample) {
  # Used to estimate sMAPE
  out_sample <- as.numeric(out_sample)
  fcast <- as.numeric(fcast)
  smape <- (abs(out_sample - fcast) * 200) / (abs(out_sample) + abs(fcast))
  return(smape)
}

mase <- function(fcast, in_sample, out_sample) {
  # Used to estimate MASE
  freq <- frequency(in_sample)
  fcast_naive_sd <- rep(NA, freq)
  for (j in (freq + 1):length(in_sample)) {
    fcast_naive_sd <- c(fcast_naive_sd, in_sample[j - freq])
  }
  masep <- mean(abs(in_sample - fcast_naive_sd), na.rm = TRUE)
  if (masep == 0) {
    print("masep == 0, setting masep = 1e-6")
    browser()
    return(1e-6)
  }

  out_sample <-
    as.numeric(out_sample)
  fcast <- as.numeric(fcast)
  mase <- (abs(out_sample - fcast)) / masep
  return(mase)
}

seasonal_naive <- function(input, fcast_horiz) {
  freq <- frequency(input)
  fcast <- naive(input, h = fcast_horiz)$mean
  if (freq > 1) {
    fcast <-
      head(rep(as.numeric(tail(input, freq)), fcast_horiz), fcast_horiz) + fcast - fcast
  }
  return(fcast)
}

theta_classic <- function(input, fcast_horiz) {
  # Used to estimate Theta classic

  # Set parameters
  wses <- wlrl <- 0.5
  theta <- 2

  # Estimate theta line (0)
  observations <- length(input)
  xt <- c(1:observations)
  xf <- c((observations + 1):(observations + fcast_horiz))
  train <- data.frame(input = input, xt = xt)
  test <- data.frame(xt = xf)

  estimate <- lm(input ~ poly(xt, 1, raw = TRUE))
  thetaline0In <- as.numeric(predict(estimate))
  thetaline0Out <- as.numeric(predict(estimate, test))

  # Estimate theta line (2)
  thetalineT <- theta * input + (1 - theta) * thetaline0In
  sesmodel <- ses(thetalineT, h = fcast_horiz)
  thetaline2In <- sesmodel$fitted
  thetaline2Out <- sesmodel$mean

  # Theta fcast
  fcastIn <- (thetaline2In * wses) + (thetaline0In * wlrl)
  fcastOut <- (thetaline2Out * wses) + (thetaline0Out * wlrl)

  # Zero fcast become positive
  for (i in 1:length(fcastOut)) {
    if (fcastOut[i] < 0) {
      fcastOut[i] <- 0
    }
  }

  return(
    list(
      fitted = fcastIn,
      mean = fcastOut,
      fitted0 = thetaline0In,
      mean0 = thetaline0Out,
      fitted2 = thetaline2In,
      mean2 = thetaline2Out
    )
  )
}

seasonality_test <- function(input, freq) {
  # Used to determine whether a time series is seasonal
  tcrit <- 1.645
  if (length(input) < 3 * freq) {
    is_seasonal <- FALSE
  } else {
    xacf <- acf(input, plot = FALSE)$acf[-1, 1, 1]
    clim <-
      tcrit / sqrt(length(input)) * sqrt(cumsum(c(1, 2 * xacf ^ 2)))
    is_seasonal <- (abs(xacf[freq]) > clim[freq])

    if (is.na(is_seasonal) == TRUE) {
      is_seasonal <- FALSE
    }
  }

  return(is_seasonal)
}

deseasonalise <- function(input, fcast_horiz) {
  # Estimate seasonally adjusted time series
  freq <- frequency(input)
  is_seasonal <- FALSE
  if (freq > 1) {
    is_seasonal <- seasonality_test(input, freq)
  }
  if (is_seasonal) {
    dec <- decompose(input, type = "multiplicative")
    output <- input / dec$seasonal
    si_out <-
      head(rep(dec$seasonal[(length(dec$seasonal) - freq + 1):length(dec$seasonal)], fcast_horiz), fcast_horiz)
  } else {
    output <- input
    si_out <- rep(1, fcast_horiz)
  }
  return(list("output" = output, "si_out" = si_out))
}

# m4_benchmarks <- function(input, fcast_horiz) {
#   # Estimate the statistical benchmarks for the M4 competition
#
#   deseason <- deseasonalise(input, fcast_horiz)
#
#   fcast_naive <- naive(input, h = fcast_horiz)$mean
#   fcast_seasonal_naive <-
#     seasonal_naive(input, fcast_horiz = fcast_horiz)
#   fcast_naive2 <-
#     naive(deseason$output, h = fcast_horiz)$mean * deseason$si_out
#   fcast_ses <-
#     ses(deseason$output, h = fcast_horiz)$mean * deseason$si_out
#   fcast_holt <-
#     holt(deseason$output, h = fcast_horiz, damped = F)$mean * deseason$si_out
#   fcast_holt_damped <-
#     holt(deseason$output, h = fcast_horiz, damped = T)$mean * deseason$si_out
#   fcast_theta_classic <-
#     theta_classic(input = deseason$output, fcast_horiz = fcast_horiz)$mean * deseason$si_out
#   fcast_combined <- (fcast_ses + fcast_holt + fcast_holt_damped) / 3
#
#   return(
#     list(
#       fcast_naive = fcast_naive,
#       fcast_seasonal_naive = fcast_seasonal_naive,
#       fcast_naive2 = fcast_naive2,
#       fcast_ses = fcast_ses,
#       fcast_holt = fcast_holt,
#       fcast_holt_damped = fcast_holt_damped,
#       fcast_theta_classic = fcast_theta_classic,
#       fcast_combined = fcast_combined
#     )
#   )
# }
