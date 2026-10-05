#' Score sleep behaviour from immobility
#'
#' This function first uses a motion classifier to decide whether an animal is moving during a given time window.
#' Then, it defines sleep as contiguous immobility for a minimum duration.
#'
#' @param data  [data.table] containing behavioural variable from or one multiple animals.
#' When it has a key, unique values, are assumed to represent unique individuals (e.g. in a [behavr] table).
#' Otherwise, it analysis the data as coming from a single animal. `data` must have a column `t` representing time.
#' @param time_window_length number of seconds to be used by the motion classifier.
#' This corresponds to the sampling period of the output data.
#' @param min_time_immobile Minimal duration (in s) of a sleep bout.
#' Immobility bouts longer or equal to this value are considered as sleep.
#' @param motion_detector_FUN function used to classify movement
#' @param untracked how windows with no tracked data enter sleep scoring. `"immobile"` (default) counts them
#' as immobility, so they can extend or create sleep bouts; `"break"` ends a bout at them, so sleep is only
#' scored where the animal was seen still. `moving` is unaffected either way. Ignored with `rule = "k"`.
#' @param rule `"classic"` (default) scores a window as sleep when no frame passed the movement threshold
#' for `min_time_immobile`. `"k"` (tentative) scores it from walking and sustained movement events,
#' ignoring tracking noise (see [sleep_rules]). It adds the columns `walking`, `sustained` and
#' `micro_awake`, leaves the classic ones as they are, never scores a window without frames as sleep,
#' and drops frames the tracker inferred. It needs 10-s windows and does not take `velocity_threshold`.
#' @param k with `rule = "k"`, sustained events within a centred 60-s window that make a window awake.
#' @param pixel with `rule = "k"`, one pixel in the units of `x`/`y`. `NULL` infers it ([pixel_size]):
#' 1 for positions in pixels, 1/500 for positions as a fraction of the ROI width (as scopr loads them).
#' @param ... extra arguments to be passed to `motion_classifier_FUN`.
#' @return a [behavr] table similar to `data` with additional variables/annotations (i.e. `moving` and `asleep`).
#' The resulting data will only have one data point every `time_window_length` seconds.
#' @details
#' The default `time_window_length` is 300 seconds -- it is also known as the "5-minute rule".
#' `sleep_annotation` is typically used for ethoscope data, whilst `sleep_dam_annotation` only works on DAM2 data.
#' These functions are *rarely used directly*, but rather passed as an argument to a data loading function,
#' so that analysis can be performed on the go.
#' @examples
# # We start by making toy data for one animal:
#' dt_one_animal <- toy_ethoscope_data(seed=2)
#' ####### Ethoscope, corrected velocity classification #########
#' sleep_dt <-  sleep_annotation(dt_one_animal, masking_duration=0)
#' print(sleep_dt)
#' # We could make a sleep `barecode'
#' \dontrun{
#' library(ggplot2)
#' ggplot(sleep_dt, aes(t,y="Animal 1",fill=asleep)) +
#'                                    geom_tile() + scale_x_time()
#' }
#' ####### Ethoscope, virutal beam cross classification #########
#' sleep_dt2 <-  sleep_annotation(dt_one_animal,
#'                              motion_detector_FUN=virtual_beam_cross_detector)
#' \dontrun{
#' library(ggplot2)
#' ggplot(sleep_dt2, aes(t,y="Animal 1",fill=asleep)) +
#'                                    geom_tile() + scale_x_time()
#' }
#' ####### DAM data, de facto beam cross classification ######
#' dt_one_animal <- toy_dam_data(seed=7)
#' sleep_dt <- sleep_dam_annotation(dt_one_animal)
#' \dontrun{
#' library(ggplot2)
#' ggplot(sleep_dt, aes(t,y="Animal 1",fill=asleep)) +
#'                                    geom_tile() + scale_x_time()
#' }
#' @seealso
#' * [motion_detectors] -- options for the `motion_detector_FUN` argument
#' * [bout_analysis] -- to further analyse sleep bouts in terms of onset and length
#' @references
#' * The relevant [rethomic tutorial section](https://rethomics.github.io/sleepr) -- on sleep analysis
#' @export
sleep_annotation <- function(data,
                            time_window_length = 10, #s
                            min_time_immobile = 300, #s = 5min
                            motion_detector_FUN = max_velocity_detector,
                            untracked = c("immobile", "break"),
                            rule = c("classic", "k"),
                            k = 3,
                            pixel = NULL,
                            ...
){
  moving = .N = is_interpolated  = .SD = asleep = NULL
  untracked <- match.arg(untracked)
  rule <- match.arg(rule)
  if(rule == "k")
    check_k_rule_args(time_window_length, k, pixel, list(...))
  # all columns likely to be needed.
  columns_to_keep <- c("t", "x", "y", "max_velocity", "velocity_threshold", "interactions",
                       "beam_crosses", "moving","asleep", "is_interpolated",
                       "walking", "sustained", "micro_awake")

  wrapped <- function(d){
    if(rule == "k")
      return(k_rule_annotation(d, time_window_length, min_time_immobile,
                               motion_detector_FUN, k, pixel, columns_to_keep, ...))
    if(nrow(d) < 100)
      return(NULL)
    # todo if t not unique, stop

    d_small <- motion_detector_FUN(d, time_window_length,...)

    if(key(d_small) != "t")
      stop("Key in output of motion_classifier_FUN MUST be `t'")

    if(nrow(d_small) < 1)
      return(NULL)
    # the times to  be queried
    time_map <- data.table::data.table(t = seq(from=d_small[1,t], to=d_small[.N,t], by=time_window_length),
                          key = "t")
    missing_val <- time_map[!d_small]

    d_small <- d_small[time_map,roll=T]
    d_small[,is_interpolated := FALSE]
    d_small[missing_val,is_interpolated:=TRUE]
    d_small[is_interpolated == T, moving := FALSE]
    sleep_breaking <- d_small$moving
    if(untracked == "break")
      sleep_breaking <- sleep_breaking | d_small$is_interpolated
    d_small[,asleep := sleep_contiguous(sleep_breaking,
                                        1/time_window_length,
                                        min_valid_time = min_time_immobile)]
    d_small <- stats::na.omit(d[d_small,
         on=c("t"),
         roll=T])
    d_small[, intersect(columns_to_keep, colnames(d_small)), with=FALSE]
  }

  if(is.null(key(data)))
     return(wrapped(data))
  data[,
       wrapped(.SD),
       by=key(data)]
}

attr(sleep_annotation, "needed_columns") <- function(motion_detector_FUN = max_velocity_detector,
                                                     rule = "classic",
                                                     ...){
  needed_columns <- attr(motion_detector_FUN, "needed_columns")
  columns <- NULL
  if(!is.null(needed_columns))
    columns <- needed_columns(...)
  # the k-rule also needs the y position of every frame
  if(identical(rule, "k"))
    columns <- unique(c(columns, "x", "y", "xy_dist_log10x1000"))
  columns
}

#' @export
#' @rdname sleep_annotation
sleep_dam_annotation <- function(data,
                                 min_time_immobile = 300){

  asleep = moving = activity = duration = .SD = . = NULL
  wrapped <- function(d){
    if(! all(c("activity", "t") %in% names(d)))
      stop("data from DAM should have a column named `activity` and one named `t`")

      out <- data.table::copy(d)
      col_order <- c(colnames(d),"moving", "asleep")
      out[, moving := activity > 0]
      bdt <- bout_analysis(moving, out)
      bdt[, asleep := duration >= min_time_immobile & !moving]
      out <- bdt[,.(t, asleep)][out, on = "t", roll=TRUE]
      data.table::setcolorder(out, col_order)
      out
  }

  if(is.null(key(data)))
    return(wrapped(data))
  data[,
       wrapped(.SD),
       by=key(data)]
}
