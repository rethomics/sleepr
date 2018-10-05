#remove data points when the time series is too sparse
curate_sparse_roi_data <- function(
  data,
  window=60,#s
  min_points=20#
){
  t_w = n_points =.N = NULL
  data[, t_w := window * floor(t/window)]
  data[,n_points := .N, by=t_w]
  data <- data[n_points > min_points]
  data[, t_w := NULL]
  data[, n_points := NULL]
  data
}
