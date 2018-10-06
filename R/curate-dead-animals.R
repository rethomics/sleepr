#' Remove -- irrelevant -- data after individual have died
#'
#' This function detects when individuals have died based on their first (very) long bout of immobility.
#' Following data (which may include spurious artefact of movement) are removed.
#'
#' @inheritParams sleep_annotation
#' @param moving_var logical variable in `data` used to define the moving (alive) state (default is `moving`)
#' @param time_window window during which to define death (default is one day)
#' @param prop_immobile proportion of immobility that counts as "dead" during time_window (see details)
#' @param resolution how much scanning windows overlap. Expressed as a factor (see details).
#'
#' @details
#' This function scans the time series looking for "death" in the right (future) data, within `time_window`.
#' Death is defined as `mean(moving_var) < prop_immobile` within a time window.
#' Moving window start every `time_window/resolution`. `resolution = 1` is fast but means no overlap.
#' The default would score an animal as dead it does not move more than *one percent of the time* for at least *one day*.
#' All data following a "death" event are removed.
#'
#' @return an object of the same type as `data` (i.e. [data.table::data.table] or [behavr::behavr]).
#'
#' @examples
#' dt1 <- toy_activity_data()
#' #all movement after day 3 is set at 0
#' dt1[t > days(3), moving := FALSE]
#' # one artefact of movement is detected at day 3.5
#' dt1[t == days(3.5), moving := TRUE]
#'
#' dt2 <- curate_dead_animals(dt1)
#' dt3 <- curate_dead_animals(dt1,prop_immobile = 0)
#' \dontrun{
#' library(ggplot2)
#' ggplot(data=dt1[,test:=1],aes(t, as.numeric(moving))) +
#'   geom_line(data=dt1[,test:=1]) +
#'   geom_line(data=dt2[, test:=2])+
#'   geom_line(data=dt3[, test:=3])+
#'   facet_grid(test ~ .)+
#'   scale_x_time()
#' }
#' @seealso
#' * [sleep_annotation] -- to score movement and slepe
# TODO
# @references
# * The relevant [rethomic tutorial section](https://rethomics.github.io/survival) -- on high-resolution survival analysis
#' @export
curate_dead_animals <- function(data,
                                moving_var = moving,
                                time_window = hours(24),
                                prop_immobile = 0.01,
                                resolution = 24){
  .N = . = moving__ = moving = .SD  = NULL
  var_name <- deparse(substitute(moving_var))
  if(!var_name %in% colnames(data))
    stop("moving_var must be a column of data. ",
         sprintf("No column named '%s'", var_name))

  wrapped <- function(do){
    d <- data.table::copy(do[, c("t", var_name),with=F])
    data.table::setnames(d, var_name, "moving__")
    target_t <-  seq(from=d[,min(t)], to= d[,max(t)], by=time_window / resolution)
    local_means <- sapply(target_t, function(x__){
      mean(d[t %between% c(x__, x__ + time_window), moving__])
    }
    )

    first_death_point <- which(local_means <= prop_immobile)[1]

    if(is.na(first_death_point))
      return(do)

    alive_dt <- data.table::data.table(t=target_t[first_death_point])
    last_valid_point <- d[moving__ == T & t < target_t[first_death_point]][.N, t]
    do[t <= last_valid_point]
  }
  if(is.null(key(data)))
    return(wrapped(data))
  data[,
       wrapped(.SD),
       by=key(data)]
}

