# devtools::install_github("ykang/tsgeneration")
library(tidyverse)
library(M4comp2018)
library(dtwclust)
# library(parallel)
library(clusterCrit)
library(ggplot2)
library(forecast)
set.seed(42)
source("m4_benchmarks_eval.R")

# TODO:
# Deseasonalise using STL

#######################################################################
# Config

m4_season <- "Monthly"
fcast_horiz <- 18
freq <- 12

num_ts <- 1000
ts_len <- 480
nrep <- 11
k_range = c(2:11)

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
    return(deseasonalise(x, fcast_horiz)))
m4_data_x_inter <-
  lapply(m4_data_x_deseason, function(ts)
    return(reinterpolate(ts$output, ts_len)))

m4_data_xx <-
  lapply(m4_data, function(ts)
    return(ts$xx))
# m4_data_type <-
#   as.integer(unlist(lapply(m4_data, function(ts)
#     return(ts$type))))
remove(m4_data)
gc(full = TRUE)

# #######################################################################
# # TS clustering
#
# # cl <- interactive_clustering(m4_data_x_deseason)
# cl_k_nrep <- tsclust(
#   m4_data_x_inter,
#   k = k_range,
#   distance = "l2",
#   centroid = "pam",
#   seed = 42,
#   trace = TRUE,
#   control = partitional_control(nrep = nrep),
#   parallel = TRUE
# )
#
# # Extract clustering results
# cl_k_nrep_k <- lapply(cl_k_nrep, function(cl)
#   return(cl@k))
# cl_k_nrep_dists <- lapply(cl_k_nrep, function(cl)
#   return(cl@cldist))
# cl_k_nrep_clusters <- lapply(cl_k_nrep, function(cl)
#   return(cl@cluster))
# remove(cl_k_nrep)
# gc(full = TRUE)
# # save(cl_k_nrep_k,
# #      cl_k_nrep_dists,
# #      cl_k_nrep_clusters,
# #      file = paste0(fname, ".RData"))
#
# #######################################################################
# # External and internal clustering metrics
#
# compute_int_metrics <- function(x) {
#   metrics <-
#     c("calinski_harabasz",
#       "gamma",
#       "gdi42",
#       "pbm",
#       "point_biserial",
#       "silhouette")
#   return(intCriteria(cl_k_nrep_dists[[x]], cl_k_nrep_clusters[[x]], metrics))
# }
#
# compute_ext_metrics <- function(x) {
#   metrics <- c("precision", "recall")
#   return(extCriteria(cl_k_nrep_clusters[[x]], m4_data_type, metrics))
# }
#
# # metrics_k_nrep <-
# #   mclapply(
# #     1:length(cl_k_nrep_k),
# #     compute_int_metrics,
# #     mc.preschedule = FALSE,
# #     mc.cores = 2,
# #     affinity.list = rep(c(2, 3), length(cl_k_nrep_k) / 2)
# #   )
# metrics_k_nrep <- lapply(1:length(cl_k_nrep_k), compute_int_metrics)
#
# metrics_df <- bind_rows(metrics_k_nrep)
# metrics_df$k <- unlist(cl_k_nrep_k)
# metrics_df %>%
#   gather(metric, value, -k) ->
#   results_df
#
# #######################################################################
# # Plot clustering metrics
#
# gg <-
#   ggplot(results_df, aes(x = k, y = value)) +
#   ggtitle(title) +
#   geom_point(size = 0.25, alpha = 0.5) +
#   geom_smooth() +
#   facet_wrap(~ metric, scales = "free")
# print(gg)
#
# ggsave(
#   paste0(fname, ".png"),
#   dpi = 100,
#   scale = 5,
#   width = 2,
#   height = 2,
#   units = "in"
# )
#
# #######################################################################
# # Extract best clusters for each k in k_range
#
# # Find median metric for each k
# metrics_df %>%
#   group_by(k) %>%
#   summarise(silhouette = median(silhouette)) ->
#   medians_df
#
# # Find row index for the median metric for each k
# metrics_df$row_id <- 1:length(cl_k_nrep_k)
# metrics_df %>%
#   select(k, silhouette, row_id) %>%
#   inner_join(medians_df) %>%
#   group_by(k) %>%
#   summarise(row = max(row_id)) ->
#   idx_df

#######################################################################
# Forecast all TS using all benchmark methods

multi_fit_ts <- function(idx) {
  fcast_naive <- naive(m4_data_x[[idx]], h = fcast_horiz)$mean
  fcast_seasonal_naive <-
    seasonal_naive(m4_data_x[[idx]], fcast_horiz = fcast_horiz)
  fcast_naive2 <-
    naive(m4_data_x_deseason[[idx]]$output, h = fcast_horiz)$mean * m4_data_x_deseason[[idx]]$si_out
  fcast_ses <-
    ses(m4_data_x_deseason[[idx]]$output, h = fcast_horiz)$mean * m4_data_x_deseason[[idx]]$si_out
  fcast_holt <-
    holt(m4_data_x_deseason[[idx]]$output,
         h = fcast_horiz,
         damped = FALSE)$mean * m4_data_x_deseason[[idx]]$si_out
  fcast_holt_damped <-
    holt(m4_data_x_deseason[[idx]]$output,
         h = fcast_horiz,
         damped = TRUE)$mean * m4_data_x_deseason[[idx]]$si_out
  fcast_theta_classic <-
    theta_classic(input = m4_data_x_deseason[[idx]]$output, fcast_horiz = fcast_horiz)$mean * m4_data_x_deseason[[idx]]$si_out

  return(
    list(
      fcast_naive = fcast_naive,
      fcast_seasonal_naive = fcast_seasonal_naive,
      fcast_naive2 = fcast_naive2,
      fcast_ses = fcast_ses,
      fcast_holt = fcast_holt,
      fcast_holt_damped = fcast_holt_damped,
      fcast_theta_classic = fcast_theta_classic,
      fcast_combined = (fcast_ses + fcast_holt + fcast_holt_damped) / 3
    )
  )
}
fcasts_m4 <- lapply(1:length(m4_data_x), multi_fit_ts)

fcast_errs <- function(idx) {
  mean_smapes <- lapply(fcasts_m4[[idx]], function(fcast) return(mean(smape(fcast, m4_data_xx[[idx]]))))
  mean_mases <- lapply(fcasts_m4[[idx]], function(fcast) return(mean(mase(fcast, m4_data_x[[idx]], m4_data_xx[[idx]]))))
  return(list(mean_smapes = mean_smapes, mean_mases = mean_mases))
}
fcast_errs_m4 <- lapply(1:length(fcasts_m4), fcast_errs)

mean_smapes <- function(fcast_name) {
  mean(unlist(lapply(1:length(fcasts_m4), function(idx) return(fcast_errs_m4[[idx]]$mean_smapes[[fcast_name]]))))
}
mean_mases <- function(fcast_name) {
  mean(unlist(lapply(1:length(fcasts_m4), function(idx) return(fcast_errs_m4[[idx]]$mean_mases[[fcast_name]]))))
}

fcast_names <- names(fcasts_m4[[1]])
mean_smapes_m4 <- lapply(fcast_names, mean_smapes)
names(mean_smapes_m4) <- fcast_names
mean_mases_m4 <- lapply(fcast_names, mean_mases)
names(mean_mases_m4) <- fcast_names

######################################################################

benchmark_names <- c("Naive", "sNaive", "Naive2", "SES", "Holt", "Damped Holt", "Theta Classic", "Combined")
smape_totals = mase_totals <- array(NA, dim = c(length(benchmark_names), fcast_horiz, length(m4_data_x)))
# Methods, Horizon, time-series
for (i in 1:length(m4_data_x)){
  in_sample <- m4_data_x[[i]]
  out_sample <- m4_data_xx[[i]]
  benchmark_fcasts <- m4_benchmarks(input=in_sample, fcast_horiz=fcast_horiz)

  # sMAPE
  for (j in 1:length(benchmark_names)){
    smape_totals[j,,i] <- smape(benchmark_fcasts[[j]], out_sample) #j the # of the benchmark
  }
  # MASE
  for (j in 1:length(benchmark_names)){
    mase_totals[j,,i] <- mase(benchmark_fcasts[[j]], in_sample, out_sample) #j the # of the benchmark
  }
}

print("########### sMAPE ###############")
for (i in 1:length(benchmark_names)){
  print(paste(benchmark_names[i], round(mean(smape_totals[i,,]), 3)))
}
print("########### MASE ################")
for (i in 1:length(benchmark_names)){
  print(paste(benchmark_names[i], round(mean(mase_totals[i,,]), 3)))
}
# print("########### OWA ################")
# for (i in 1:length(benchmark_names)){
#   print(paste(benchmark_names[i],
#               round(((mean(mase_totals[i,,])/mean(mase_totals[3,,]))+(mean(smape_totals[i,,])/mean(smape_totals[3,,])))/2, 3)))
# }

# fit_cluster_ts <- function(cl_number, cl_assignment) {
#   cl_idxes <- match(cl_assignment, cl_number)
#   print(paste0("Cluster #", cl_number, " has length: ", sum(cl_idxes), na.rm = TRUE))
#
#   fcast_cl <- lapply(cl_idxes, fit_ts)
# }
#
# get_clusters <- function(row_idx) {
#   k <- cl_k_nrep_k[[row_idx]]
#   cl_assignment <- cl_k_nrep_clusters[[row_idx]]
#   print(paste0("Number of TS per cluster for k=", k))
#   print(table(cl_assignment))
#
#   res <- lapply(1:k, fit_clusters, cl_assignment)
# }
#
# res <- lapply(idx_df$row, get_clusters)
