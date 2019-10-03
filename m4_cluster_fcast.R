# devtools::install_github("ykang/tsgeneration")
library(tidyverse)
library(M4comp2018)
library(dtwclust)
library(ggplot2)

set.seed(42)
options(warn = 2, width = 1024)
source("fcast.R")
source("cluster.R")

###########################################################################
# Config ####

m4_season <- "Monthly"
fcast_horiz <- 18
freq <- 12

num_ts <- NA #46341
ts_len <- 480
nrep <- 7
k_range <- c(2+1:20*2)
err_names <- c("sMAPE", "MASE", "OWA")

###########################################################################
# Preprocess M4 data ####

if (is.na(num_ts)) {
  m4_data <- Filter(function(ts)
    ts$period == m4_season, M4)
} else {
  m4_data <- sample(Filter(function(ts)
    ts$period == m4_season, M4), num_ts)
}

fname <-
  paste0(
    "nts",
    length(m4_data),
    "_m4_",
    tolower(substr(m4_season, 1, 3)),
    "_tslen",
    ts_len,
    "_med_l2_nrep",
    nrep
  )
title <-
  paste0(
    length(m4_data),
    " TS from M4 ",
    m4_season,
    ", interpolated to length ",
    ts_len,
    ", Median + L2, clustering from k=",
    min(k_range),
    " to k=",
    max(k_range),
    ", ",
    nrep,
    " clustering reps:"
  )
print(title)

# M4 Competition only data
m4_data_x <-
  lapply(m4_data, function(ts)
    return(ts$x))
m4_data_xx <-
  lapply(m4_data, function(ts)
    return(ts$xx))

# Deseasonalise and linearly interpolate to ts_len
m4_data_x_deseason <-
  lapply(m4_data_x, function(x)
    return(deseasonalise(x, fcast_horiz)))
m4_data_x_inter <-
  lapply(m4_data_x_deseason, function(ts)
    return(reinterpolate(ts$output, ts_len)))

remove(m4_data)
gc(verbose = TRUE)
print(summary(unlist(m4_data_x)))

###########################################################################
# Forecast each TS using each of the benchmark methods ####

print("M4 Competition data:")
fcasts <-
  lapply(1:length(m4_data_x),
         multi_fit_ts,
         m4_data_x,
         m4_data_x_deseason)
fcast_names <- names(fcasts[[1]])

###########################################################################
# Compute sMAPE, MASE, and OWA ####

fcast_errs <-
  lapply(1:length(fcasts),
         compute_fcast_errs,
         fcasts,
         m4_data_x,
         m4_data_xx)
mean_errs_df <-
  as.data.frame(lapply(fcast_names, mean_fcast_errs, fcast_errs))
colnames(mean_errs_df) <- fcast_names
mean_errs_df <-
  rbind(mean_errs_df,
        colMeans(mean_errs_df / mean_errs_df$naive2))
rownames(mean_errs_df) <- err_names
print(round(mean_errs_df, 3))
write_csv(as.data.frame(t(mean_errs_df)), path = paste0("benchmark_m4_", m4_season, ".csv"))

# ###########################################################################
# # Cluster M4 deseasonalised and interpolated TS ####
#
# print(paste0(
#   "Clustering from k=",
#   min(k_range),
#   " to k=",
#   max(k_range),
#   ", for ",
#   nrep,
#   " reps"
# ))
# cl <- cluster_ts(m4_data_x_inter, k_range, nrep)
#
# print(paste0(
#   "Computing clustering metrics for ",
#   length(cl$k_nrep_k),
#   " cluster models"
# ))
# cl_metrics_df <- compute_cl_metrics(cl)
#
# ###########################################################################
# # Select the best forecast type per M4 TS cluster ####
#
# print("Finding best clustered forecasts based on OWA:")
# # Find median metric for each k
# cl_metrics_df %>%
#   group_by(k) %>%
#   summarise(pbm = median(pbm)) ->
#   cl_medians_df
#
# # Find row index for the median metric for each k
# cl_metrics_df$row <- 1:length(cl$k_nrep_k)
# cl_metrics_df %>%
#   select(k, pbm, row) %>%
#   inner_join(cl_medians_df) %>%
#   group_by(k) %>%
#   summarise(row = max(row)) ->
#   idx_df
#
# cl_best_ks <-
#   lapply(
#     idx_df$row,
#     find_best_clusters,
#     cl,
#     fcast_names,
#     fcast_errs,
#     mean_errs_df$naive2[1:2]
#   )
# names(cl_best_ks) <- idx_df$k
# # print(lapply(cl_best_ks, round, 3))
#
# ###########################################################################
# # Plot clustering errors as a function of k_range ####
#
# cl_best_ks %>%
#   bind_rows %>%
#   t %>%
#   as.data.frame ->
#   cl_best_ks_df
# cl_best_ks_df$k <- idx_df$k
# colnames(cl_best_ks_df) <- c(err_names, "k")
# print("Best OWA clustered result:")
# print(round(cl_best_ks_df[which.min(cl_best_ks_df$OWA), ], 3))
#
# cl_best_ks_df %>%
#   gather(metric, error,-k) ->
#   results_df
#
# benchmark_best_df <-
#   data.frame(metric = err_names,
#              error = unlist(lapply(err_names, function(err)
#                return(
#                  min(mean_errs_df[err, ])
#                ))))
#
# gg <-
#   ggplot(results_df, aes(x = k, y = error)) +
#   ggtitle(title) +
#   geom_line() +
#   geom_point() +
#   facet_wrap( ~ metric, scales = "free_y") +
#   geom_hline(data = benchmark_best_df, aes(yintercept = error), colour = "red")
# print(gg)
# ggsave(
#   paste0("stmod_", fname, ".png"),
#   dpi = 100,
#   scale = 5,
#   width = 2,
#   height = 2,
#   units = "in"
# )
