#' Find "bouts" in categorical time series
#'
#' This function is used to find contiguous regions of unique value
#' in a -- potentially irregular/heterogeneous -- univariate categorical time series.
#'
#' @param var name of the variable to use from `data`
#' @inheritParams sleep_annotation
#' @return an object of the same type as `data` (i.e. [data.table::data.table] or [behavr::behavr]).
#' Each row is a specific bout characterised by three columns.
#' * `t` -- its *onset*
#' * `duration` --  its length
#' * `<var>` -- a column with the same name as `var`. The value of `var` for this bout.
#' @examples
#' #TODO
#' @seealso
#' todo
# * Tutorial on anlysis of sleep architecture \url{http://gilestrolab.github.io/rethomics/tutorial/sleep_quantification.html#sleep-architecture} TODO
# * [rle] run length encoding function, on which this analysis is based
#' @export
bout_analysis <- function(var,data){

  var_name <- deparse(substitute(var))
  wrapped <- function(d){
    dt <- data.table::copy(d[, c("t",var_name), with=F])
    data.table::setnames(dt,var_name,"var__")

    dt[,delta_t:= c(0,diff(dt[,t]))]
    r <- rle(dt[, var__])
    vals <-r$values
    r$values <- 1:length(r$values)

    bdt <- data.table::data.table(delta_t = dt[,delta_t],
                                  bout_id__ = inverse.rle(r),
                                  key="bout_id__")

    bout_times <- bdt[,
                      .(duration = sum(delta_t)),
                      by="bout_id__"]
    r$values <- vals #?
    dt <- data.table::data.table(
      var__ = vals,
      duration = bout_times[,duration],
      t = cumsum(c(0,bout_times[1:.N-1,duration])) + dt[1,t]
    )
    data.table::setnames(dt,"var__",var_name)
    dt
  }

  if(is.null(key(data)))
    return(wrapped(data))
  data[,
       wrapped(.SD),
       by=key(data)]
}

#
# make_bout_dt <- function(x, sdt){
#   #sdt <- copy(sub_data)
#   sdt[,delta_t:= c(0,diff(sdt[,t]))]
#   r <- rle(x)
#   vals <-r$values
#   r$values <- 1:length(r$values)
#
#   bdt <- data.table::data.table(delta_t = sdt[,delta_t], ..bout_id..=inverse.rle(r),key="..bout_id..")
#
#   bout_times <- bdt[,
#                     .(duration = sum(delta_t)),
#                     by="..bout_id.."]
#   r$values <- vals
#   out <- data.table::data.table(
#     x = vals,
#     duration = bout_times[,duration],
#     t = cumsum(c(0,bout_times[1:.N-1,duration])) + sdt[1,t]
#   )
#   var_name <- deparse(substitute(var))
#   setnames(out,"x",var_name)
#   out
# }
