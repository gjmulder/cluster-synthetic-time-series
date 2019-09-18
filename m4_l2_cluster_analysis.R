# devtools::install_github("ykang/tsgeneration")
library(tidyr)
library(dplyr)
library(purrr)
library(M4comp2018)
library(dtwclust)
library(parallel)
library(clusterCrit)
library(ggplot2)

set.seed(42)

# TODO:
# Deseasonalise using STL

#######################################################################
# Config

num_ts <- 20000
ts_len <- 240
reps <- 10
k_range = c(2:19)

title <-
  paste0(
    num_ts,
    " TS sampled from M4 monthly, reinterpolated to ",
    ts_len,
    " length, ",
    reps,
    " k-means iteration(s)"
  )
fname <-
  paste0("nts", num_ts, "_m4-mon_tslen", ts_len, "_reps", reps)

data(M4)
monthly_m4 <-
  sample(Filter(function(l)
    l$period == "Monthly", M4), num_ts)
print(summary(unlist(lapply(monthly_m4, function(x)
  return(x$n)))))
monthly_m4_inter <-
  lapply(monthly_m4, function(ts)
    return(reinterpolate(ts$x, ts_len)))
remove(M4, monthly_m4)
gc(full = TRUE)

#######################################################################
# TS clustering

# cl <- interactive_clustering(tsl2)

# Cluster
cl_k_reps <- tsclust(
  monthly_m4_inter,
  # type = "partitional",
  k = k_range,
  # preproc = "NULL",
  distance = "L2",
  centroid = "pam",
  seed = 42,
  trace = TRUE,
  control = partitional_control(nrep = reps),
  parallel = TRUE
)

# Extract clustering results
cl_k_reps_k <- lapply(cl_k_reps, function(cl)
  return(cl@k))
cl_k_reps_dists <- lapply(cl_k_reps, function(cl)
  return(cl@cldist))
cl_k_reps_clusters <- lapply(cl_k_reps, function(cl)
  return(cl@cluster))
remove(cl_k_reps)
gc(full = TRUE)

#######################################################################
# Internal clustering metrics
compute_metrics <- function(x) {
  metrics <-
    c(
      "calinski_harabasz",
      "davies_bouldin",
      "gamma",
      "gdi42",
      "pbm",
      "point_biserial",
      "sd_dis",
      "silhouette",
      "wemmert_gancarski"
    )
  return(intCriteria(cl_k_reps_dists[[x]], cl_k_reps_clusters[[x]], "all"))
}
metrics_k_reps <-
  mclapply(
    1:length(cl_k_reps_k),
    compute_metrics,
    mc.preschedule = FALSE,
    mc.cores = 2,
    affinity.list = rep(c(2, 3), length(cl_k_reps_k) / 2)
  )
# metrics_k_reps <- lapply(1:length(cl_k_reps_k), compute_metrics)

int_metrics_df <- bind_rows(metrics_k_reps)
int_metrics_df$k <- unlist(cl_k_reps_k)
int_metrics_df %>%
  gather(metric, value, -k) %>%
  filter(metric %in% metrics) ->
  int_results_df

#######################################################################
# Plot clustering metrics

gg <-
  ggplot(int_results_df, aes(x = k, y = value)) +
  ggtitle(title) +
  geom_point(size = 0.25, alpha = 0.5) +
  geom_smooth() +
  facet_wrap(~ metric, scales = "free")
print(gg)

ggsave(
  paste0(fname, ".png"),
  dpi = 100,
  scale = 4,
  width = 2,
  height = 2,
  units = "in"
)
