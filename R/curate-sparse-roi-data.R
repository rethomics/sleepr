#remove data points when the time series is too sparse
curate_sparse_roi_data <- function(
  data,
  window=60,#s
  min_points=20#
){

  d <- copy(data)
  d[, t_w := window * floor(t/window)]
  sparsity <- d[, t_w := window * floor(t/window)]
  d[,sparsity := .N,by=t_w]
  d <- d[sparsity >min_points,]
  d$t_w <- NULL
  d$sparsity <- NULL
  d
}
