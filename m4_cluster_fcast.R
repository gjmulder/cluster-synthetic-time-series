# devtools::install_github("ykang/tsgeneration")
library(tidyverse)
library(M4comp2018)
library(dtwclust)
library(ggplot2)

set.seed(42)
options(warn = 2)
source("fcast.R")
source("cluster.R")

###########################################################################
# Config ####

m4_season <- "Monthly"
fcast_horiz <- 18
freq <- 12

num_ts <- NA
ts_len <- 480
nrep <- 11
k_range <- c(3:20)
err_names <- c("sMAPE", "MASE", "OWA")

title <-
  paste0(
    num_ts,
    " TS sampled from M4 deseasonalised ",
    m4_season,
    ", ",
    ts_len,
    " interpolated length, PAM + L2, clustering from k=",
    min(k_range),
    " to k=",
    max(k_range),
    ", ",
    nrep,
    " clustering reps:"
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

###########################################################################
# Preprocess M4 data ####

print(title)
m4_data <-
  Filter(function(ts)
    ts$period == m4_season, M4)
# sample(Filter(function(ts)
#   ts$period == m4_season, M4), num_ts)

# M4 Competition only data
m4_data_x <-
  lapply(m4_data, function(ts)
    return(subset(ts$x, end = length(ts$x) - fcast_horiz)))
m4_data_xx <-
  lapply(m4_data, function(ts)
    return(subset(ts$x, start = (
      length(ts$x) - fcast_horiz + 1
    ))))

# Deseasonalise and linearly interpolate to ts_len
m4_data_x_deseason <-
  lapply(m4_data_x, function(x)
    return(deseasonalise(x, fcast_horiz)))
m4_data_x_inter <-
  lapply(m4_data_x_deseason, function(ts)
    return(reinterpolate(ts$output, ts_len)))

# m4_data_type <-
#   as.integer(unlist(lapply(m4_data, function(ts)
#     return(ts$type))))

print(summary(unlist(lapply(m4_data, function(x)
  return(x$n)))))

###########################################################################
# Cluster M4 deseasonalised and interpolated TS ####

print(paste0(
  "Clustering from k=",
  min(k_range),
  " to k=",
  max(k_range),
  ", for ",
  nrep,
  " reps"
))
cl <- cluster_ts(m4_data_x_inter, k_range, nrep)

print(paste0(
  "Computing clustering metrics for ",
  length(cl$k_nrep_k),
  " cluster models"
))
cl_metrics_df <- compute_cl_metrics(cl)

###########################################################################
# Post-M4 Competition data ####

m4_data_post_x <-
  lapply(m4_data, function(ts)
    return(ts$x))
m4_data_post_xx <-
  lapply(m4_data, function(ts)
    return(ts$xx))
m4_data_x_post_deseason <-
  lapply(m4_data_post_x, function(x)
    return(deseasonalise(x, fcast_horiz)))

remove(m4_data)
gc(verbose = TRUE)

###########################################################################
# Forecast each TS using each of the benchmark methods ####

print("M4 Competition estimates:")
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

###########################################################################
# Post-M4 forecast each TS using each of the benchmark methods ####

print("M4 Competition benchmark results:")
fcasts_post <-
  lapply(1:length(m4_data_post_x),
         multi_fit_ts,
         m4_data_post_x,
         m4_data_x_post_deseason)

###########################################################################
# Post M4 compute sMAPE, MASE, and OWA ####

fcast_post_errs <-
  lapply(
    1:length(fcasts_post),
    compute_fcast_errs,
    fcasts_post,
    m4_data_post_x,
    m4_data_post_xx
  )
mean_errs_post_df <-
  as.data.frame(lapply(fcast_names, mean_fcast_errs, fcast_post_errs))
colnames(mean_errs_post_df) <- fcast_names
mean_errs_post_df <-
  rbind(mean_errs_post_df,
        colMeans(mean_errs_post_df / mean_errs_post_df$naive2))
rownames(mean_errs_post_df) <- err_names
print(round(mean_errs_post_df, 3))

###########################################################################
# Select the best forecast type per M4 TS cluster ####

cl_select_best_fcast <-
  function(cl_n,
           cl_assignment,
           fcast_names,
           fcast_post_errs) {
    cl_assignment %>% when(cl_assignment == cl_n) -> cl_n_match
    print(paste0("Cluster #", cl_n, " has size: ", sum(cl_n_match)))
    cl_n_idx <- c(1:length(cl_assignment))[cl_n_match]

    # Compute mean errors for cl_n
    cl_fcast_errs_df <-
      as.data.frame(lapply(fcast_names, mean_fcast_errs, fcast_errs[cl_n_idx]))
    colnames(cl_fcast_errs_df) <- fcast_names
    cl_fcast_errs_df <-
      rbind(cl_fcast_errs_df,
            colMeans(cl_fcast_errs_df / cl_fcast_errs_df$naive2))
    rownames(cl_fcast_errs_df) <- err_names
    # print(cl_fcast_errs_df)

    # Find best forecast method using OWA for cl_n
    best_cl_n <- names(which.min(cl_fcast_errs_df["OWA",]))
    print(best_cl_n)

    # Return the out of sample errors
    best_cl_err <-
      bind_cols(lapply(fcast_post_errs[cl_n_idx], function(fcast_post_err, best_err)
        return(fcast_post_err[best_err]), best_cl_n))
    return(best_cl_err)
  }

get_clusters <-
  function(idx,
           cl,
           fcast_names,
           fcast_post_errs,
           naive2_errs_post) {
    k <- cl$k_nrep_k[[idx]]
    cl_assignment <- cl$k_nrep_clusters[[idx]]
    print("")
    print(paste0("k=", k))
    cl_best <-
      rowMeans(bind_cols(
        lapply(
          1:k,
          cl_select_best_fcast,
          cl_assignment,
          fcast_names,
          fcast_post_errs
        )
      ))
    cl_best_v <- c(cl_best, mean(cl_best / naive2_errs_post))
    names(cl_best_v) <- err_names
    # print(cl_best_v)
    return(cl_best_v)
  }

print("Finding best clustered forecasts based on OWA:")
# Find median metric for each k
cl_metrics_df %>%
  group_by(k) %>%
  summarise(pbm = median(pbm)) ->
  cl_medians_df

# Find row index for the median metric for each k
cl_metrics_df$row <- 1:length(cl$k_nrep_k)
cl_metrics_df %>%
  select(k, pbm, row) %>%
  inner_join(cl_medians_df) %>%
  group_by(k) %>%
  summarise(row = max(row)) ->
  idx_df

cl_best_ks <-
  lapply(
    idx_df$row,
    get_clusters,
    cl,
    fcast_names,
    fcast_post_errs,
    mean_errs_post_df$naive2[1:2]
  )
names(cl_best_ks) <- idx_df$k
# print(lapply(cl_best_ks, round, 3))

###########################################################################
# Plot clustering errors as a function of k_range ####

cl_best_ks %>%
  bind_rows %>%
  t %>%
  as.data.frame ->
  cl_best_ks_df
cl_best_ks_df$k <- idx_df$k
colnames(cl_best_ks_df) <- c(err_names, "k")
print("Best OWA clustered result:")
print(round(cl_best_ks_df[which.min(cl_best_ks_df$OWA),], 3))

cl_best_ks_df %>%
  gather(metric, error, -k) ->
  results_df

benchmark_best_df <-
  data.frame(metric = err_names,
             error = unlist(lapply(err_names, function(err)
               return(
                 min(mean_errs_post_df[err,])
               ))))

gg <-
  ggplot(results_df, aes(x = k, y = error)) +
  ggtitle(title) +
  geom_line() +
  facet_wrap(~ metric, scales = "free_y") +
  geom_hline(data = benchmark_best_df, aes(yintercept = error), colour = "red")
print(gg)
ggsave(
  paste0(fname, ".png"),
  dpi = 100,
  scale = 5,
  width = 2,
  height = 2,
  units = "in"
)
