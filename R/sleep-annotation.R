#' @importFrom data.table ":="
#' @importFrom data.table "key"
NULL

#' Infer sleep from immobility bouts
#'
#' This function first uses a motion classifier to decide whether an animal is moving during a given time window.
#' Then, it defines sleep as contiguous immobility for a minimal duration.
#'
#' @param data  [behavr] table containing behavioural variable from multiple (or one) animals.
#' It must contain, at least, the column `t` in addition to the variables needed for motion classification.
#' @param time_window_length number of seconds to be used by the motion classifier.
#' This corresponds to the sampling period of the output data.
#' @param min_time_immobile Minimal duration (in s) of a sleep bout.
#' Immobility bouts longer or equal to this value are considered as sleep.
#' @param motion_detector_FUN function used to classify movement
#' @param ... extra arguments to be passed to `motion_classifier_FUN`.
#' @return a [behavr] table similar to `data` with additional variables/annotations (i.e. `moving` and `asleep`).
#' The resulting data will only have one data point every `time_window_length` seconds.
#' @details
#' The default `time_window_length` is 300 seconds also known as the "5min rule".
#' `sleepAnnotation` is typically used for ethoscope data, whilst `sleepDAMAnnotation` only works on DAM data.
#' These function is *rarely used directly*, but rather passed as an argument to a data loading function,
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
# @seealso
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
  d <- copy(data)
  ori_keys <- key(d)
  d <- curateSparseRoiData(d)
  if(nrow(d) <1)
    return(NULL)

  d_small <- motion_classifier_FUN(d, time_window_length,...)

  # special variable "has interacted". We sum over it
  if("has_interacted" %in% colnames(d)){
    d_n_interations <- d[, .(n_interactions = sum(has_interacted)), by=key(d)]
    d_small <- d_small[d_n_interations]
    d[, has_interacted := NULL]
  }

  if(key(d_small) != "t_round")
    stop("Key in output of motion_classifier_FUN MUST be `t_round'")
  setnames(d_small,"t_round", "t")
  d$t <- NULL
  d_small <- d_small[unique(d,by=key(d))]

  t_out <- seq(from=d_small[1,t], to=d_small[.N,t], by=time_window_length)

  time_map <- data.table(t=t_out,key="t")
  missing_val <- time_map[!d_small]

  d_small <- d_small[time_map,roll=T]
  d_small[,is_interpolated := FALSE]
  d_small[missing_val,is_interpolated:=TRUE]

  d_small[is_interpolated == T, moving := FALSE]
  d_small[,asleep := sleep_contiguous(moving,
                                      1/time_window_length,
                                      min_valid_time = min_time_immobile)]
  setkeyv(d_small, ori_keys)
  na.omit(d_small)
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

  moving_dt[, asleep := sleepr:::sleep_contiguous(moving, 1/60, 300), by=key(data)]
  moving_dt[data, on=c(key(data), "t")]
}
