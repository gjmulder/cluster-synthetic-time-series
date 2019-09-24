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

fcast_horiz <- 18
freq <- 12 #The frequency of the data

num_ts <- 1000
ts_len <- 500
nrep <- 5
k_range = c(5:8)

title <-
  paste0(num_ts,
         " TS sampled from M4, PAM + L2, ",
         ts_len,
         " length, ",
         nrep,
         " clustering reps")
fname <-
  paste0("nts",
         num_ts,
         "_m4_pam_l2_intmet_tslen_fitcl",
         ts_len,
         "_nrep",
         nrep)

m4_data <-
  sample(Filter(function(ts)
    ts$period == "Monthly", M4), num_ts)
# sample(M4, num_ts)
print(summary(unlist(lapply(m4_data, function(x)
  return(x$n)))))
m4_data_x <-
  lapply(m4_data, function(ts)
    return(ts$x))
m4_data_x_deseason <-
  lapply(m4_data_x, deseasonalise, fcast_horiz)

# m4_data_x_inter <-
#   lapply(m4_data, function(ts)
#   return(reinterpolate(ts$x, ts_len)))
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

# cl <- interactive_clustering(tsl2)

# Cluster
cl_k_nrep <- tsclust(
  m4_data_x_deseason,
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
  facet_wrap( ~ metric, scales = "free")
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
#
# cl_fcasts <- function(input, fh) {
#   #Used to estimate the statistical benchmarks of the M4 competition
#
#   res <- deseasonalise(input, fh)
#   des_input <- res[[1]]
#   SIout <- res[[2]]
#
#   f1 <- naive(input, h = fh)$mean #Naive
#   f2 <- naive_seasonal(input, fh = fh) #Seasonal Naive
#   f3 <- naive(des_input, h = fh)$mean * SIout #Naive2
#   f4 <- ses(des_input, h = fh)$mean * SIout #Ses
#   f5 <- holt(des_input, h = fh, damped = F)$mean * SIout #Holt
#   f6 <- holt(des_input, h = fh, damped = T)$mean * SIout #Damped
#   f7 <-
#     Theta.classic(input = des_input, fh = fh)$mean * SIout #Theta
#   f8 <- (f4 + f5 + f6) / 3 #Comb
#
#   return(list(f1, f2, f3, f4, f5, f6, f7, f8))
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
