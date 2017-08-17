#' Finds "bouts" in categorical time series.
#'
#' This function is used to find contiguous regions of unique value
#' in a -- potentially irregular -- univariate categorical time series.
#'
#' @param var Column variable to use in \code{data}
#' @param data  `data.table` representing a behaviour of multiple animals.
#' @return A `data.table` where **each row is a bout** and with the columns:
#' * `start_time` The onset of this bout (second, same reference as `t`)
#' * `length` The bout duration (second)
#' * A column with the same name as `var` that describe the value of var for this whole bout.
#' @details
#' Bout analysis will be performed by biological individual.
#' Individuals are defined by the `data.table` key.
#' @examples
#' set.seed(1)
#' # 1000 points the first 500 points should have higher chance to be 1 than the last 500:
#' y_var <- round(c(runif(500,0,1),
#'                    runif(500,0,0.75)))
#' # first 500 point are for individual "A", next 500 points are for "B":
#' dt <- data.table( y = y_var,
#'              t = rep(1:500,2)*12,
#'              id = rep(c("A","B"),each=500),key="id")
#' print(dt)
#' bout_dt <- boutAnalysis(y,dt)
#' print(bout_dt)
#' # average duration and number of bouts for both individuals and
#' # bouts of y=0 or y=1
#' summary_dt <- bout_dt[,
#'         .(n=.N,
#'           mean_duration=mean(length))
#'         ,by=c(key(bout_dt),"y")]
#' print(summary_dt)
#' @seealso
#' * Tutorial on anlysis of sleep architecture \url{http://gilestrolab.github.io/rethomics/tutorial/sleep_quantification.html#sleep-architecture} TODO
#' * [rle] run length encoding function, on which this analysis is based
#' @export
bout_analysis <- function(var,data){
  dt <- copy(as.data.table(data))
  var_name <- deparse(substitute(var))
  setnames(dt,var_name,"var")

  out <- dt[,makeBoutDt(var,.SD),by=c(key(dt))]

  setnames(out,"var",var_name)
  out
}


makeBoutDt <- function(x,sub_data){
  sdt <- copy(sub_data)
  sdt[,delta_t:= c(0,diff(sub_data[,t]))]
  r <- rle(x)
  vals <-r$values
  r$values <- 1:length(r$values)

  bdt <- data.table(delta_t = sdt[,delta_t], bout_id=inverse.rle(r),key="bout_id")

  bout_times <- bdt[,list(length=sum(delta_t)),by="bout_id"]
  r$values <- vals
  out <- data.table(
    x = vals,
    length = bout_times[,length],
    start_time = cumsum(c(0,bout_times[1:.N-1,length])) + sdt[1,t]
  )
  var_name <- deparse(substitute(var))
  setnames(out,"x",var_name)
  out
}
