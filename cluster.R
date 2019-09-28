library(dtwclust)
library(clusterCrit)
library(ggplot2)
library(parallel)

#######################################################################
# TS clustering

# cl <- interactive_clustering(m4_data_x_deseason)
cluster_ts <- function(data_x, k_range, nrep) {
  cl_k_nrep <- tsclust(
    data_x,
    k = k_range,
    distance = "l2",
    centroid = "pam",
    seed = 42,
    trace = FALSE,
    control = partitional_control(nrep = nrep),
    parallel = TRUE
  )

  # Extract clustering results
  k_nrep_k <- lapply(cl_k_nrep, function(cl)
    return(cl@k))
  k_nrep_dists <- lapply(cl_k_nrep, function(cl)
    return(cl@cldist))
  k_nrep_clusters <- lapply(cl_k_nrep, function(cl)
    return(cl@cluster))

  remove(cl_k_nrep)
  gc(full = TRUE)

  # save(cl_k_nrep_k,
  #      cl_k_nrep_dists,
  #      cl_k_nrep_clusters,
  #      file = paste0(fname, ".RData"))
  return(
    list(
      k_nrep_k = k_nrep_k,
      k_nrep_dists = k_nrep_dists,
      k_nrep_clusters = k_nrep_clusters
    )
  )
}

#######################################################################
# External and internal clustering metrics

compute_cl_metrics <-
  function(cl) {
    compute_int_metrics <- function(x) {
      metrics <- c("pbm")
        # c("calinski_harabasz",
        #   "gamma",
        #   "gdi42",
        #   "pbm",
        #   "point_biserial",
        #   "silhouette")
      return(intCriteria(cl$k_nrep_dists[[x]], cl$k_nrep_clusters[[x]], metrics))
    }

    # compute_ext_metrics <- function(x) {
    #   metrics <- c("precision", "recall")
    #   return(extCriteria(cl$k_nrep_clusters[[x]], m4_data_type, metrics))
    # }

    # metrics_k_nrep <-
    #   mclapply(
    #     1:length(cl$k_nrep_k),
    #     compute_int_metrics,
    #     mc.preschedule = FALSE,
    #     mc.cores = 2,
    #     affinity.list = rep(c(2, 3), length(cl$k_nrep_k) / 2)
    #   )
    metrics_k_nrep <-
      lapply(1:length(cl$k_nrep_k), compute_int_metrics)

    metrics_df <- bind_rows(metrics_k_nrep)
    metrics_df$k <- unlist(cl$k_nrep_k)

    return(metrics_df)
  }

plot_metrics <- function(metrics_df, title, fname) {
  metrics_df %>%
    gather(metric, value, -k) ->
    results_df

  gg <-
    ggplot(results_df, aes(x = k, y = value)) +
    ggtitle(title) +
    geom_point(size = 0.25, alpha = 0.5) +
    geom_smooth() +
    facet_wrap(~ metric, scales = "free")
  print(gg)

  ggsave(
    paste0(fname, ".png"),
    dpi = 100,
    scale = 5,
    width = 2,
    height = 2,
    units = "in"
  )
}
