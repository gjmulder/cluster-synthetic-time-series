# devtools::install_github("ykang/tsgeneration")
# library(doParallel)
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
# Add stochiastic noise

#######################################################################
#
# Config variables

num_ts <- 10
# vals <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29)
# selected_features <- c('entropy', 'trend', 'seasonal_strength')
# extra_target <- c(0.6, 0.2, 0.1)
vals_df <- data.frame("entropy" = c(1:3) / 4, "trend" = c(1:3) / 4)
targets_df <- as.data.frame(t(expand.grid(vals_df)))
extra_targets <- NULL
# ground_truths <- unlist(lapply(1:length(vals), rep, num_ts))

#######################################################################

fname <-
  paste0(
    "numts=",
    num_ts * nrow(targets_df),
    "_num_feat=",
    nrow(targets_df),
    "_feat=",
    paste0(selected_features, collapse = ","),
    "tgt=",
    paste0(extra_targets, collapse = ",")
  )

#######################################################################
# Generate synthetic time series

# cl <- makeCluster(2)
# registerDoParallel(cl)

gen_ts <- function(targets, n) {
  print(paste0("val=(",
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
      features = c('entropy', 'stl_features'),
      selected.features = colnames(vals_df),
      target = c(targets, extra_target)
    )
  )
  names(df_ts) <-
    paste0("t", paste0(targets, collapse = "-"), "_", c(1:n))
  return(df_ts)
}

tsl <- lapply(vals_comb_df, gen_ts, num_ts)
# # stopCluster(cl)
# # registerDoSEQ()

ts_df <- do.call("cbind", tsl)
tsl2 <- as.list(ts_df)

save(tsl2, file = paste0(fname, ".Rdata"))

#######################################################################
# DTW clustering

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
k_range = c(2:(2 * length(vals)))
metrics_k <- lapply(k_range, ts_clust_k)
# stopCluster(cl)
# registerDoSEQ()

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
    num_ts * nrow(targets_df),
    " TS, name=",
    colnames(targets_df),
    ", num feats=",
    nrow(targets_df),
    paste0(colnames(vals_df), collapse = ","),
    "), target=(",
    paste0(extra_target, collapse = ","),
    ")"
  )

gg <-
  ggplot(int_results) +
  ggtitle(title) +
  geom_point(aes(x = k, y = value), size = 0.5) +
  geom_vline(xintercept = length(vals)) +
  facet_wrap(~ metric, scales = "free")

print(gg)
ggsave(
  paste0(fname, ".png"),
  scale = 4,
  width = 10,
  height = 10,
  units = "cm"
)
