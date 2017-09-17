#' Score sleep behaviour from immobility
#'
#' This function first uses a motion classifier to decide whether an animal is moving during a given time window.
#' Then, it defines sleep as contiguous immobility for a minimal duration.
#'
#' @param data  [data.table] containing behavioural variable from or one multiple animals.
#' When it has a key, unique values, are assumed to represent unique inviduals (e.g. in a [behavr] table).
#' Otherwise, it analysis the data as comming from a single animal. `data` must have a column `t` representing time.
#' @param time_window_length number of seconds to be used by the motion classifier.
#' This corresponds to the sampling period of the output data.
#' @param min_time_immobile Minimal duration (in s) of a sleep bout.
#' Immobility bouts longer or equal to this value are considered as sleep.
#' @param motion_detector_FUN function used to classify movement
#' @param ... extra arguments to be passed to `motion_classifier_FUN`.
#' @return a [behavr] table similar to `data` with additional variables/annotations (i.e. `moving` and `asleep`).
#' The resulting data will only have one data point every `time_window_length` seconds.
#' @details
#'
#' The default `time_window_length` is 300 seconds also known as the "5 minute rule".
#' `sleep_annotation` is typically used for ethoscope data, whilst `sleep_dam_annotation` only works on DAM2 data.
#' These functions are *rarely used directly*, but rather passed as an argument to a data loading function,
#' so that analysis can be performed on the go.
#' @examples
#' #todo
# # We strat by making toy data for one animal:
# dt_one_animal <- toyEthoscopeData(seed=2)
# ####### Ethoscope, corrected velocity classification #########
# sleep_dt <-  sleepAnnotation(dt_one_animal, masking_duration=0)
# print(sleep_dt)
# # We make a sleep `barecode'
# ggplot(sleep_dt, aes(t,y="Animal 1",fill=asleep)) +
#                                    geom_tile() + scale_x_time()
# ####### Ethoscope, virutal beam cross classification #########
# sleep_dt2 <-  sleepAnnotation(dt_one_animal,
#                              motion_classifier_FUN=virtualBeamCrossClassif)
# ggplot(sleep_dt, aes(t,y="Animal 1",fill=asleep)) +
#                                    geom_tile() + scale_x_time()
# #' ####### DAM data, de facto beam cross classification ######
# dt_one_animal <- toyDAMData(seed=7)
# sleep_dt <- sleepDAMAnnotation(dt_one_animal)
# ggplot(sleep_dt, aes(t,y="Animal 1",fill=asleep)) +
#                                    geom_tile() + scale_x_time()
#' @seealso
#' * [motion_detectors] -- options for the `motion_detector_FUN` argument
#' * [bout_analysis] -- to further analyse sleep bouts in terms of onset and length
# * Tutorial for sleep analysis with ethoscopes \url{http://gilestrolab.github.io/rethomics/tutorial/todo}
# * [loadEthoscopeData] and [loadDAM2Data] to load data and optionally apply these analysis on the fly.
# * [maxVelocityClassifierMasked] the motion classifiers that can be used.
#' @export
sleep_annotation <- function(data,
                            time_window_length = 10, #s
                            min_time_immobile = 300, #s = 5min
                            motion_detector_FUN = max_velocity_detector,
                            ...
){
  # d <- copy(data)
  # ori_keys <- key(d)
  # d <- curateSparseRoiData(d)
  # todo warn?

  # all columns likely to be needed.
  columns_to_keep <- c("t", "x", "y", "max_velocity", "interactions",
                       "beam_crosses", "moving","asleep", "is_interpolated")

  wrapped <- function(d){
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
    d_small[,asleep := sleep_contiguous(moving,
                                        1/time_window_length,
                                        min_valid_time = min_time_immobile)]
    d_small <- na.omit(d[d_small,
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
                                                     ...){
  needed_columns <- attr(motion_detector_FUN, "needed_columns")
  if(!is.null(needed_columns))
    needed_columns(...)
}

#' @export
#' @rdname sleep_annotation
sleep_dam_annotation <- function(
  data,
  time_window_length = 60, #s
  min_time_immobile = 300 # s 5min rule
){
  # if(!is.null(motion_classifier_FUN))
  #   stop("Cannot use a motion classifier with DAM data")
  if(! behavr::is.behavr(data))
    stop("data should be a `behavr` table!")

  if(! all(c("activity", "t") %in% names(data)))
    stop("data from DAM should have a column named `activity` and one named `t`")

  # todo here, check for irregular time series
  # and concistancy between diff(t) and time_window_length


  moving_dt <- data[,
                     .(t=t,
                       moving = activity>0),
                     by=key(data)
                    ]

  moving_dt[, asleep := sleep_contiguous(moving, 1/60, 300), by=key(data)]
  moving_dt[data, on=c(key(data), "t")]
}
