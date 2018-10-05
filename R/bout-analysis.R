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
#' # Bout analysis on binary variable:
#' dt <- toy_dam_data()
#' dt[, moving := activity > 0]
#' bdt <- bout_analysis(moving,dt)
#' print(bdt)
#' # With multiple states
#' dt <- toy_ethoscope_data()
#' # we discretise x position in three states: left, middle and right (1/3 each)
#' dt[, location := as.character( cut(x,
#'                                breaks = c(0.0, .33, .67, 1.0),
#'                                labels = c("left", "middle", "right")))]
#'
#' bdt <- bout_analysis(location, dt)
#' @seealso
#' * [sleep_annotation] -- to generate a binary sleep variable
#' * [rle] run length encoding function -- on which this analysis is based
#' @references
#' * The relevant [rethomic tutorial section](https://rethomics.github.io/sleepr) -- on sleep analysis
#' @export
bout_analysis <- function(var,data){
  .SD = NULL
  var_name <- deparse(substitute(var))
  if(!var_name %in% colnames(data))
    stop("var must be a column of data. ",
          sprintf("No column named '%s'", var_name))
  if(is.null(key(data)))
    return(boot_analysis_wrapped(data, var_name))
  data[,
       boot_analysis_wrapped(.SD, var_name),
       by=key(data)]
}

#' @noRd
boot_analysis_wrapped <- function(d, var_name){
  var__ = . = delta_t = bout_id__ = duration = .N  = NULL
  dt <- data.table::copy(d[, c("t",var_name), with=F])
  data.table::setnames(dt,var_name,"var__")
  dt[,delta_t:= c(diff(dt[,t]), 0)]
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
