# devtools::install_github("ykang/tsgeneration")
library(tidyr)
library(dplyr)
library(purrr)
library(M4comp2018)
library(dtwclust)
library(parallel)
library(clusterCrit)
library(ggplot2)
library(viridis)

set.seed(42)

# TODO:
# Add single seasonality and deseasonalise using STL

#######################################################################
# Config variables

ts_len <- 240
iters <- 3
desc <-
  paste0("M4 monthly dataset, reinterpolated to ",
         ts_len,
         " length, ",
         iters,
         " iteration(s)")
short_desc <- paste0("m4-mon_tslen", ts_len, "_iters", iters)

fname <- short_desc
k_range = c(2:19)

# data("M4")
monthly_m4 <-
  sample(Filter(function(l)
    l$period == "Monthly", M4), 10000)
print(summary(unlist(lapply(monthly_m4, function(x)
  return(x$n)))))
monthly_m4_inter <-
  lapply(monthly_m4, function(ts)
    return(reinterpolate(ts$x, ts_len)))

#######################################################################
# TS clustering

# cl <- interactive_clustering(tsl2)

# Cluster and validate
cl_k <- tsclust(
  monthly_m4_inter,
  # type = "partitional",
  k = k_range,
  # preproc = "NULL",
  distance = "L2",
  centroid = "pam",
  seed = 42,
  trace = TRUE,
  control = partitional_control(nrep = iters),
  parallel = TRUE
)
gc(full = TRUE)

# Internal clustering metrics, i.e. no ground truths
compute_metrics <- function(cl) {
  metrics <-
    c(
      "ball_hall",
      "banfeld_raftery",
      "c_index",
      "calinski_harabasz",
      "davies_bouldin",
      "det_ratio",
      "pbm",
      "point_biserial",
      "ratkowsky_lance",
      # "ray_turi",
      "sd_dis",
      "sd_scat",
      "silhouette"
    )
  return(intCriteria(cl@cldist, cl@cluster, metrics))
}
# metrics_k <-
#   mclapply(
#     cl_k,
#     compute_metrics,
#     mc.preschedule = FALSE,
#     mc.cores = 2,
#     affinity.list = rep(c(2, 3), length(cl_k) / 2)
#   )

cl_idx <- 1:(iters * length(k_range))
cl_dists <- lapply(cl_idx, function(x) return(cl_k[[x]]@cldist))
cl_clusters <- lapply(cl_idx, function(x) return(cl_k[[x]]@cluster))
# iter_idx <- unlist(lapply(c(1:iters), rep, length(k_range)))

int_metrics_df <- bind_rows(metrics_k)
int_metrics_df$k <- rep(k_range, iters)
int_metrics_df$iter <- unlist(lapply(c(1:iters), rep, length(k_range)))
int_metrics_df %>%
  gather(metric, value, -k, -iter) %>%
  mutate(iter = as.character(iter)) ->
  # filter(iter == 1) %>%
  # arrange(iter, k) ->
  int_results

#######################################################################
# Plot clustering metrics

title <-
  paste0("Num TS=",
         length(monthly_m4_inter),
         " ",
         desc)

gg <-
  ggplot(int_results) +
  ggtitle(title) +
  geom_line(aes(x = k, y = value, colour = iter), size = 0.5) +
  # scale_color_viridis(discrete = TRUE) +
  # geom_vline(xintercept = ncol(targets_df)) +
  facet_wrap( ~ metric, scales = "free")
print(gg)

ggsave(
  paste0(paste0(
    "nts=",
    length(monthly_m4_inter),
    "_",
    fname,
    ".png"
  )),
  dpi = 100,
  scale = 5,
  width = 3,
  height = 2,
  units = "in"
)
