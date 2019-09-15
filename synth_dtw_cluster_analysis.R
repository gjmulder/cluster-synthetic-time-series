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

desc <- "reduce feature range by 10X"
short_desc <- "reduc-feat"

num_ts <- 10
# selected_features <- c('entropy', 'trend', 'seasonal_strength')
# extra_target <- c(0.6, 0.2, 0.1)
noise <- 0.4
features <- c('entropy', 'stl_features')
# target_ranges_df <-
#   data.frame(
#     "linearity"   = c(0, 5, 10),
#     "trend"       = c(0:2) / 3,
#     "spike"       = c(0:2) / 3,
#     "curvature"   = c(-5, 0, 5)
#   )

target_ranges_df <-
  data.frame(
    "linearity"   = c(0, 5, 10) / 10.0,
    "trend"       = 0.5 + c(-1:1) / 10.0,
    "spike"       = 0.5 + c(-1:1) / 10.0,
    "curvature"   = c(-5, 0, 5) / 10.0
  )
targets_df <- as.data.frame(t(expand.grid(target_ranges_df)))
targets_df <-
  rbind(targets_df, sample(trunc(num_ts / 2):trunc(3 * num_ts / 2), size =
                             ncol(targets_df), replace = TRUE))

extra_features <- c("entropy")
extra_targets <- c(noise)

# ground_truths <- unlist(lapply(1:length(vals), rep, num_ts))

#######################################################################
# Generate synthetic time series

gen_ts <- function(targets_n) {
  targets <- head(targets_n, -1)
  n <- tail(targets_n, 1)
  print(paste0("targets=(",
               paste0(targets, collapse = ","),
               "), n=",
               n))

  df_ts <- as.data.frame(
    generate_ts_with_target(
      parallel = TRUE,
      # seed = 42,
      n = n,
      ts.length = 100,
      freq = 12,
      seasonal = 1,
      features = features,
      selected.features = c(colnames(target_ranges_df), extra_features),
      target = c(targets, extra_targets)
    )
  )
  names(df_ts) <-
    paste0("t", paste0(targets, collapse = "-"), "_", 1:n)
  return(df_ts)
}

tsl <- lapply(targets_df, gen_ts)
ts_df <- do.call("cbind", tsl)
tsl2 <- as.list(ts_df)

######################################################################

synth_data_fname <-
  paste0(
    short_desc,
    "_nfeat=",
    ncol(targets_df),
    "_lfeat=",
    nrow(target_ranges_df),
    "_feat=",
    paste0(rownames(targets_df), collapse = ","),
    "_efeat=",
    paste0(extra_features, collapse = ","),
    "_etgt=",
    paste0(extra_targets, collapse = ",")
  )

save(
  desc,
  short_desc,
  tsl2,
  target_ranges_df,
  extra_features,
  extra_targets,
  file = paste0(synth_data_fname, ".Rdata")
)
# load(paste0(synth_data_fname, ".Rdata"))

# add_noise_zscore <- function(ts, noise) {
#   magnitude <- max(abs(ts))
#   noise_ts <- runif(length(ts),
#                     min = (noise * -magnitude),
#                     max = (noise * magnitude))
#   # return(zscore(ts + noise_ts))
#   return(ts + noise_ts)
# }
# ts_noise <- lapply(tsl2, add_noise_zscore, noise)

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
    tsl2,
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
    "Num TS=",
    length(tsl2),
    ", uniq feats=",
    ncol(targets_df),
    ", feat len=",
    nrow(target_ranges_df),
    ", tgts=(",
    paste0(rownames(targets_df), collapse = ","),
    "), +feat=(",
    paste0(extra_features, collapse = ","),
    "), +tgts=(",
    paste0(extra_targets, collapse = ","),
    ")"
  )

gg <-
  ggplot(int_results) +
  ggtitle(title) +
  geom_point(aes(x = k, y = value), size = 0.5) +
  geom_vline(xintercept = ncol(targets_df)) +
  facet_wrap( ~ metric, scales = "free")
print(gg)

ggsave(
  paste0(paste0(
    "nts=",
    length(tsl2),
    "_",
    synth_data_fname,
    ".png"
  )),
  scale = 4,
  width = 10,
  height = 10,
  units = "cm"
)
