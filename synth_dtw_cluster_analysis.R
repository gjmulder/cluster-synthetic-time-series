# devtools::install_github("ykang/tsgeneration")
library(doParallel)
library(tidyr)
library(dplyr)
library(tsgeneration)
library(dtwclust)
library(clusterCrit)
library(ggplot2)

set.seed(42)

num_ts <- 50
# freqs <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29)
# selected_features <- c('entropy', 'trend', 'seasonal_strength')
# target <- c(0.6, 0.2, 0.1)
freqs <- c(2:11)
selected_features <- c('entropy', 'seasonal_strength')
target <- c(0.6, 0.001)

fname <-
  paste0(
    "num",
    num_ts * length(freqs),
    "_freq=",
    paste0(freqs, collapse = ","),
    "_feat=",
    paste0(selected_features, collapse = ","),
    "_tgt=",
    paste0(target, collapse = ",")
  )
#######################################################################
# Generate synthetic time series

cl <- makeCluster(2)
registerDoParallel(cl)

gen_ts <- function(freq, n) {
  print(paste0("freq=", freq, ", n=", n))
  df_freq <- as.data.frame(
    generate_ts_with_target(
      parallel = TRUE,
      # seed = 42,
      n = n,
      ts.length = 100,
      freq = freq,
      seasonal = 1,
      features = c('entropy', 'stl_features'),
      selected.features = selected_features,
      target = target
    )
  )
  names(df_freq) <- paste0("f", freq, "_", c(1:n))
  return(df_freq)
}

tsl <- lapply(freqs, gen_ts, num_ts)
# # stopCluster(cl)
# # registerDoSEQ()

ts_df <- do.call("cbind", tsl)
tsl2 <- as.list(ts_df)

save(tsl2, file = paste0(fname, ".Rdata"))

#######################################################################
# DTW clustering

# cl <- makeCluster(2)
# registerDoParallel(cl)

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
  return(cl@cluster)
}

k_range = c(2:(2 * length(freqs)))
cl_k_df <- as.data.frame(sapply(k_range, ts_clust_k))
names(cl_k_df) <- k_range

# stopCluster(cl)
# registerDoSEQ()

#######################################################################
# Validate clustering qulaity

val_criteria <- function(vec)
  extCriteria(vec, ground_truths, "all")

ground_truths <- unlist(lapply(1:length(freqs), rep, num_ts))
ext_metrics_df <-
  as.data.frame(do.call(rbind, lapply(cl_k_df, val_criteria)))
ext_metrics_df$k <- as.numeric(rownames(ext_metrics_df))

ext_metrics_df %>%
  gather(metric, value, -k) %>%
  mutate(value = unlist(value)) ->
  results

#######################################################################
# Plot clustering metrics

title <-
  paste0(
    num_ts * length(freqs),
    " TS, freqs=(",
    paste0(freqs, collapse = ","),
    "), feats=(",
    paste0(selected_features, collapse = ","),
    "), tgt=(",
    paste0(target, collapse = ","),
    ")"
  )

gg <-
  ggplot(results) +
  ggtitle(title) +
  geom_point(aes(x = k, y = value)) +
  geom_vline(xintercept = length(freqs)) +
  facet_wrap( ~ metric, scales = "free")

print(gg)
ggsave(
  paste0(fname, ".png"),
  scale = 2,
  width = 20,
  height = 10,
  units = "cm"
)
