# devtools::install_github("ykang/tsgeneration")
library(tidyr)
library(dplyr)
library(purrr)
library(tsgeneration)
library(dtwclust)
library(clusterCrit)
library(ggplot2)

set.seed(42)

# TODO:
# Add single seasonality and deseasonalise using STL

#######################################################################
# Config variables

num_ts <- 50
# selected_features <- c('entropy', 'trend', 'seasonal_strength')
# extra_target <- c(0.6, 0.2, 0.1)
noise <- 0.4
features <- c('entropy', 'stl_features')
vals_comb_df <-
  data.frame(
    "entropy" = c(1:4) / 5,
    "trend" = c(0:3) / 4,
    "spike" = c(0:3) / 4
  )
targets_df <- as.data.frame(t(expand.grid(vals_comb_df)))
extra_targets <- NULL
# ground_truths <- unlist(lapply(1:length(vals), rep, num_ts))

#######################################################################
# Generate synthetic time series

gen_ts <- function(targets, n) {
  n_rand <- sample(trunc(n/2):trunc(3*n/2), 1)
  print(paste0("targets=(",
               paste0(targets, collapse = ","),
               "), n=",
               n_rand))

  df_ts <- as.data.frame(
    generate_ts_with_target(
      parallel = TRUE,
      # seed = 42,
      n = n_rand,
      ts.length = 100,
      freq = 12,
      seasonal = 1,
      features = features,
      selected.features = colnames(vals_comb_df),
      target = c(targets, extra_targets)
    )
  )
  names(df_ts) <-
    paste0("t", paste0(targets, collapse = "-"), "_", c(1:n))
  return(df_ts)
}

tsl <- lapply(targets_df, gen_ts, num_ts)
ts_df <- do.call("cbind", tsl)
tsl2 <- as.list(ts_df)

#######################################################################

synth_data_fname <-
  paste0(
    "nfeat=",
    ncol(targets_df),
    "_feat=",
    paste0(rownames(targets_df), collapse = ","),
    "_ext=",
    paste0(extra_targets, collapse = ","),
    ".Rdata"
  )
# save(tsl2, synth_data_fname)
load(synth_data_fname)

add_noise_zscore <- function(ts, noise) {
  magnitude <- max(abs(ts))
  noise_ts <- runif(
    length(ts),
    min = (noise * -magnitude),
    max = (noise * magnitude)
  )
  # return(zscore(ts + noise_ts))
  return(ts + noise_ts)
}
ts_noise <- lapply(tsl2, add_noise_zscore, noise)

#######################################################################
# TS clustering

# cl <- makeCluster(2)
# registerDoParallel(cl)
ext_criteria <- function(vec)
  extCriteria(vec, ground_truths, "all")

int_criteria <- function(traj, vec)
  intCriteria(traj, vec, "all")

# cl <- interactive_clustering(tsl2)
ts_clust_k <- function(k) {
  cl <- tsclust(
    ts_noise,
    k = k,
    preproc = NULL,
    distance = "L2",
    centroid = "pam",
    seed = 42,
    trace = TRUE,
    control = partitional_control(nrep = 1),
    parallel = FALSE
  )

  # Validate clustering
  ext_metrics = list() # ext_criteria(cl@cluster)
  int_metrics = int_criteria(cl@cldist, cl@cluster)
  return(list(ext_metrics, int_metrics))
}

# Cluster and Validate
k_range = c(2:(2 * length(targets_df)))
metrics_k <- lapply(k_range, ts_clust_k)

# # External metrics, i.e. with known ground truths
# ext_metrics_df <- bind_rows(lapply(metrics_k, `[[`, 1))
# ext_metrics_df$k = k_range
# ext_metrics_df %>%
#   gather(metric, value, -k) ->
#   ext_results

# Internal metrics, i.e. no ground truths
int_metrics_df <- bind_rows(lapply(metrics_k, `[[`, 2))
int_metrics_df$k = k_range
int_metrics_df %>%
  gather(metric, value, -k) ->
  int_results

#######################################################################
# Plot clustering metrics

title <-
  paste0(
    num_ts * ncol(targets_df),
    " TS, noise=",
    noise,
    ", uniq feats=",
    ncol(targets_df),
    ", targets=(",
    paste0(rownames(targets_df), collapse = ","),
    "), extra targets=(",
    paste0(extra_targets, collapse = ","),
    ")"
  )

gg <-
  ggplot(int_results) +
  ggtitle(title) +
  geom_point(aes(x = k, y = value), size = 0.5) +
  geom_vline(xintercept = ncol(targets_df)) +
  facet_wrap(~ metric, scales = "free")
print(gg)

ggplot_fname <-
  paste0(
    "nts=",
    length(tsl2),
    "_noise=",
    noise,
    "_nfeat=",
    ncol(targets_df),
    "_feat=",
    paste0(rownames(targets_df), collapse = ","),
    "_ext=",
    paste0(extra_targets, collapse = ","),
    ".png"
  )
ggsave(
  paste0(ggplot_fname),
  scale = 4,
  width = 10,
  height = 10,
  units = "cm"
)
