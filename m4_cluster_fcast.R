# devtools::install_github("ykang/tsgeneration")
library(tidyverse)
library(M4comp2018)
library(dtwclust)
# library(parallel)
library(clusterCrit)
library(ggplot2)
# library(forecast)
set.seed(42)
source("m4_benchmarks_eval.R")

# TODO:
# Deseasonalise using STL

#######################################################################
# Config

m4_season <- "Monthly"
fcast_horiz <- 18
freq <- 12

num_ts <- 5000
ts_len <- 500
nrep <- 5
k_range = c(2:19) # c(5:8)

title <-
  paste0(
    num_ts,
    " TS sampled from M4 deseasonalised ",
    m4_season,
    ", ",
    ts_len,
    " interpolated length, PAM + L2, ",
    nrep,
    " clustering reps"
  )
fname <-
  paste0(
    "nts",
    num_ts,
    "_m4_",
    tolower(substr(m4_season, 1, 3)),
    "_tslen",
    ts_len,
    "_pam_l2_intmet_nrep",
    nrep
  )

#######################################################################
# Preprocess M4 data

m4_data <-
  sample(Filter(function(ts)
    ts$period == m4_season, M4), num_ts)
# sample(M4, num_ts)
print(summary(unlist(lapply(m4_data, function(x)
  return(x$n)))))

m4_data_x <-
  lapply(m4_data, function(ts)
    return(ts$x))
m4_data_x_deseason <-
  lapply(m4_data_x, function(x)
    return(deseasonalise(x, fcast_horiz)$des_input))
m4_data_x_inter <-
  lapply(m4_data_x_deseason, function(ts)
  return(reinterpolate(ts, ts_len)))

m4_data_xx <-
  lapply(m4_data, function(ts)
    return(ts$xx))
# m4_data_type <-
#   as.integer(unlist(lapply(m4_data, function(ts)
#     return(ts$type))))
remove(m4_data)
gc(full = TRUE)

#######################################################################
# TS clustering

# cl <- interactive_clustering(m4_data_x_deseason)
cl_k_nrep <- tsclust(
  m4_data_x_inter,
  k = k_range,
  distance = "l2",
  centroid = "pam",
  seed = 42,
  trace = TRUE,
  control = partitional_control(nrep = nrep),
  parallel = TRUE
)

# Extract clustering results
cl_k_nrep_k <- lapply(cl_k_nrep, function(cl)
  return(cl@k))
cl_k_nrep_dists <- lapply(cl_k_nrep, function(cl)
  return(cl@cldist))
cl_k_nrep_clusters <- lapply(cl_k_nrep, function(cl)
  return(cl@cluster))
remove(cl_k_nrep)
gc(full = TRUE)
save(cl_k_nrep_k,
     cl_k_nrep_dists,
     cl_k_nrep_clusters,
     file = paste0(fname, ".RData"))

#######################################################################
# External and internal clustering metrics

compute_int_metrics <- function(x) {
  metrics <-
    c(
      "calinski_harabasz",
      "gamma",
      "gdi42",
      "pbm",
      "point_biserial",
      "silhouette",
      "wemmert_gancarski"
    )
  return(intCriteria(cl_k_nrep_dists[[x]], cl_k_nrep_clusters[[x]], metrics))
}

compute_ext_metrics <- function(x) {
  metrics <- c("precision", "recall")
  return(extCriteria(cl_k_nrep_clusters[[x]], m4_data_type, metrics))
}

# metrics_k_nrep <-
#   mclapply(
#     1:length(cl_k_nrep_k),
#     compute_int_metrics,
#     mc.preschedule = FALSE,
#     mc.cores = 2,
#     affinity.list = rep(c(2, 3), length(cl_k_nrep_k) / 2)
#   )
metrics_k_nrep <- lapply(1:length(cl_k_nrep_k), compute_int_metrics)

metrics_df <- bind_rows(metrics_k_nrep)
metrics_df$k <- unlist(cl_k_nrep_k)
metrics_df %>%
  gather(metric, value, -k) ->
  results_df

#######################################################################
# Plot clustering metrics

gg <-
  ggplot(results_df, aes(x = k, y = value)) +
  ggtitle(title) +
  geom_point(size = 0.25, alpha = 0.5) +
  geom_smooth() +
  facet_wrap(~ metric, scales = "free")
print(gg)

ggsave(
  paste0(fname, ".png"),
  dpi = 100,
  scale = 5,
  width = 2,
  height = 2,
  units = "in"
)

# ######################################################################
# m4_benchmarks <- function(input, fcast_horiz) {
#   # Estimate the statistical benchmarks for the M4 competition
#
#   deseason <- deseasonalise(input, fcast_horiz)
#
#   fcast_naive <- naive(input, h = fcast_horiz)$mean
#   fcast_seasonal_naive <-
#     seasonal_naive(input, fcast_horiz = fcast_horiz)
#   fcast_naive2 <-
#     naive(deseason$des_input, h = fcast_horiz)$mean * deseason$si_out
#   fcast_ses <-
#     ses(deseason$des_input, h = fcast_horiz)$mean * deseason$si_out
#   fcast_holt <-
#     holt(deseason$des_input, h = fcast_horiz, damped = F)$mean * deseason$si_out
#   fcast_holt_damped <-
#     holt(deseason$des_input, h = fcast_horiz, damped = T)$mean * deseason$si_out
#   fcast_theta_classic <-
#     theta_classic(input = deseason$des_input, fcast_horiz = fcast_horiz)$mean * deseason$si_out
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
#
# ######################################################################
#
# benchmark_names <- c("Naive", "sNaive", "Naive2", "SES", "Holt", "Damped Holt", "Theta Classic", "Combined")
# smape_totals = mase_totals <- array(NA, dim = c(length(benchmark_names), fcast_horiz, length(m4_data_x)))
# # Methods, Horizon, time-series
# for (i in 1:length(m4_data_x)){
#   in_sample <- m4_data_x[[i]]
#   out_sample <- m4_data_xx[[i]]
#   benchmark_fcasts <- m4_benchmarks(input=in_sample, fcast_horiz=fcast_horiz)
#
#   # sMAPE
#   for (j in 1:length(benchmark_names)){
#     smape_totals[j,,i] <- smape_cal(out_sample, benchmark_fcasts[[j]]) #j the # of the benchmark
#   }
#   # MASE
#   for (j in 1:length(benchmark_names)){
#     mase_totals[j,,i] <- mase_cal(in_sample, out_sample, benchmark_fcasts[[j]]) #j the # of the benchmark
#   }
# }
#
# #################################################################################
#
# print("########### sMAPE ###############")
# for (i in 1:length(benchmark_names)){
#   print(paste(benchmark_names[i], round(mean(smape_totals[i,,]), 3)))
# }
# print("########### MASE ################")
# for (i in 1:length(benchmark_names)){
#   print(paste(benchmark_names[i], round(mean(mase_totals[i,,]), 3)))
# }
# print("########### OWA ################")
# for (i in 1:length(benchmark_names)){
#   print(paste(benchmark_names[i],
#               round(((mean(mase_totals[i,,])/mean(mase_totals[3,,]))+(mean(smape_totals[i,,])/mean(smape_totals[3,,])))/2, 3)))
# }
