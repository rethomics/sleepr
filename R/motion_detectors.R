#' Motion detector for Ethocope data
#'
#' Defines whether a *single animal* is moving according to:
#'
#' * Validated and corrected subpixel velocity ([max_velocity_detector]), the most rigorous
#' * Uncorrected subpixel velocity ([max_velocity_detector_legacy])
#' * Crossing a virtual beam in the middle of the region of interest ([virtual_beam_cross_detector])
#'
#' [max_velocity_detector] is the default movement classification for real-time ethoscope experiments.
#' It is benchmarked against human-generated ground truth.
#' @name motion_detectors
#' @param data [data.table::data.table] containing behavioural variables of *a single animal* (no id).
#' It must have the columns `xy_dist_log10x1000`(for computing subpixel velocity),
#' `x`(beam cross), `t` and `has_interacted` (whether a stimulus was delivered).
#' @param velocity_correction_coef an empirical coefficient to correct velocity with respect
#'  to variable framerate.
#' @inheritParams sleep_annotation
#' @param masking_duration number of second during which any movement is ignored (velocity is set to 0) after
#' a stimulus is delivered (aka interaction).
#' @param velocity_threshold uncorrected velocity above which an animal is classified as `moving' (for the legacy version).
#' @return an object of the same type as `data` (i.e. [data.table::data.table] or [behavr::behavr])  with additional columns:
#' * `moving` Logical, TRUE iff. motion was detected.
#' * `beam_crosses` The number of beam crosses
#' (when the animal crosses x = 0.5 -- that is the midpoint of the region of interest) within the time window
#' * `max_velocity` The maximal velocity within the time window.
#' The resulting data is sampled at a period equals to `time_window_length`.
#' @details
#'  These functions are *rarely called directly*, but typically used is in the context of [sleep_annotation].
#' @seealso
#' TODO
#' * [sleep_annotation] -- which requieres a motion detector
#' @export
max_velocity_detector  <- function(data,
                                   time_window_length,
                                   velocity_correction_coef =3e-3,
                                   masking_duration=6){

  d <- prepare_data_for_motion_detector(data,
                                        c("t", "xy_dist_log10x1000", "x"),
                                        time_window_length,
                                        "has_interacted")
  d[,dt := c(NA,diff(t))]
  #d[,surface_change := xor_dist * 1e-3]
  d[,dist := 10^(xy_dist_log10x1000/1000) ]
  d[,velocity := dist/dt]

  a = velocity_correction_coef

  d[,beam_cross := abs(c(0,diff(sign(.5 - x))))]
  d[,beam_cross := as.logical(beam_cross)]

  # masking here
  if(!"has_interacted" %in% colnames(d)){
    if(masking_duration >0)
      warning("Data does not contain an `has_interacted` column.
              Cannot apply masking!.
              Set `masking_duration = 0` to ignore masking")
    d[, has_interacted := 0]
  }

  d[,interaction_id := cumsum(has_interacted)]
  d[,
    masked := t < (t[1] + masking_duration),
    by=interaction_id
    ]
  d[ ,velocity := ifelse(masked & interaction_id != 0, 0, velocity)]
  d[,beam_cross := !masked & beam_cross]
  d[,interaction_id := NULL]
  d[,masked := NULL]
  # end of masking

  d[, velocity_corrected :=  velocity  * dt  /a]
  d_small <- d[,.(
    max_velocity = max(velocity_corrected[2:.N]),
    # dist = sum(dist[2:.N]),
    interactions = as.integer(sum(has_interacted)),
    beam_crosses = as.integer(sum(beam_cross))
  ), by="t_round"]

  d_small[, moving :=  ifelse(max_velocity > 1, TRUE,FALSE)]
  data.table::setnames(d_small, "t_round", "t")
  d_small
}

#' @export
#' @rdname motion_detectors
max_velocity_detector_legacy <- function(data, velocity_threshold=.006){

  d <- prepare_data_for_motion_detector(data,
                                        c("t", "xy_dist_log10x1000"),
                                        time_window_length)

  d[,dt := c(NA,diff(t))]
  d[,velocity := 10^(xy_dist_log10x1000/1000)/dt ]
  #d[,max_velocity := 10^(xy_dist_log10x1000/1000)/dt ]
  d_small <- d[,.(
    max_velocity = max(velocity)
  ), by="t_round"]

  d_small[, moving :=  ifelse(max_velocity > velocity_threshold, TRUE,FALSE)]
  data.table::setnames(d_small, "t_round", "t")
  d_small
}


#' @export
#' @rdname motion_detectors
virtual_beam_cross_detector <- function(data, time_window_length){
  d <- prepare_data_for_motion_detector(data,
                                        c("t", "x"),
                                        time_window_length)
  d[,beam_cross := abs(c(0,diff(sign(.5 - x))))]
  d[,beam_cross := as.logical(beam_cross)]

  d_small <- d[,
               .(moving = any(beam_cross)),
               by="t_round"]
  data.table::setnames(d_small, "t_round", "t")
  d_small
}

#' copy needed columns, and add t_round var for downsampling
#' @noRd
prepare_data_for_motion_detector <- function(data,
                                             needed_columns,
                                             time_window_length,
                                             optional_columns=NULL){
  # todo assert no key/unique
  if(! all(needed_columns %in% names(data)))
    stop(sprintf("data from ethoscope should have columns named %s!", paste(needed_columns, collapse=", ")))
  needed_columns <- unique(c(needed_columns, intersect(names(data),optional_columns)))
  d <- data.table::copy(data[, needed_columns, with=FALSE])
  d[, t_round := time_window_length * floor(t /time_window_length)]
  d <- curate_sparse_roi_data(d)
  data.table::setkeyv(d, "t_round")
}


