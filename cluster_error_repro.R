library(dtwclust)
for (num_ts in 46342:46341) {
  print(num_ts)
  cl_k_nrep <-
    tsclust(lapply(1:num_ts, function(x)
      return(0)),
      k = 2,
      distance = "l2")
}
