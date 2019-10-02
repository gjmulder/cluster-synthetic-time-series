library(dtwclust)
library(clusterCrit)
library(ggplot2)
library(parallel)

#######################################################################
# TS clustering

cluster_ts <- function(data_x, k_range, nrep) {
  cl_k_nrep <- tsclust(
    data_x,
    k = k_range,
    distance = "l2",
    centroid = "median",
    seed = 42,
    trace = TRUE,
    control = partitional_control(nrep = nrep,
                                  iter.max = 1000),
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
  gc(verbose = TRUE)

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

###########################################################################
# Select the best forecast type per M4 TS cluster ####

cl_select_best_fcast <-
  function(cl_n,
           cl_assignment,
           fcast_names,
           fcast_errs) {
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
    best_cl_n <- names(which.min(cl_fcast_errs_df["OWA", ]))
    print(best_cl_n)

    # Return the best errors
    best_cl_err <-
      bind_cols(lapply(fcast_errs[cl_n_idx], function(fcast_err, best_err)
        return(fcast_err[best_err]), best_cl_n))
    return(best_cl_err)
  }

find_best_clusters <-
  function(idx,
           cl,
           fcast_names,
           fcast_errs,
           naive2_errs) {
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
          fcast_errs
        )
      ))
    cl_best_v <- c(cl_best, mean(cl_best / naive2_errs))
    names(cl_best_v) <- err_names
    # print(cl_best_v)
    return(cl_best_v)
  }

###########################################################################
# Plot clustering metrics as a function of k_range ####

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
    paste0("intmet_", fname, ".png"),
    dpi = 100,
    scale = 5,
    width = 2,
    height = 2,
    units = "in"
  )
}
