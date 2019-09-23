# devtools::install_github("ykang/tsgeneration")
library(tidyverse)
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
# ts_len <- 100
nrep <- 3
k_range = c(2:19)

title <-
  paste0(
    num_ts,
    " TS sampled from M4, PAM + L2, reinterpolated to ",
    # " TS sampled from M4, TS type as ground truth, reinterpolated to ",
    # ts_len,
    # " length, ",
    nrep,
    " clustering reps"
  )
fname <-
  paste0("nts", num_ts, "_m4_pam_l2_intmet_nrep", nrep)
# paste0("nts", num_ts, "_m4_pam_intmet_tslen", ts_len, "_nrep", nrep)
m4_data <-
  sample(M4, num_ts)
# ts$period == "Monthly" && ts$type == "Finance", M4), num_ts)
print(summary(unlist(lapply(m4_data, function(x)
  return(x$n)))))
m4_data_x <-
  lapply(m4_data, function(ts)
    return(ts$x))
# return(reinterpolate(ts$x, ts_len)))
m4_data_type <-
  as.integer(unlist(lapply(m4_data, function(ts)
    return(ts$type))))
remove(m4_data)
gc(full = TRUE)

#######################################################################
# TS clustering

# cl <- interactive_clustering(tsl2)

# Cluster
cl_k_nrep <- tsclust(
  m4_data_x,
  # type = "partitional",
  k = k_range,
  # preproc = "NULL",
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
