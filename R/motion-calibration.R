#' Tracking-noise diagnostics
#'
#' [max_velocity_detector] calls a time window "moving" when its peak corrected
#' velocity exceeds a threshold. When a motionless animal's tracking noise sits
#' near it, spurious "movements" fragment every immobility bout: a bout of 5 min
#' needs 30 consecutive immobile 10-s windows, so a per-window false positive rate
#' `p` leaves only `(1 - p)^30` of true rests intact.
#'
#' These helpers support [motion_qc], which reports how often a still animal
#' crosses the fixed threshold. Windows in which the animal is still are found
#' *from position alone* -- its position over the half-minute either side of the
#' window stays within one pixel -- so the check does not depend on the velocity
#' it checks. Where the fixed threshold fails, score sleep with
#' `sleep_annotation(rule = "k")` ([sleep_rules]). The algorithms and defaults are
#' identical to ethoscopy's `motion_calibration` module.
#'
#' @name motion_calibration
#' @param d_small binned data of *a single animal*, with columns `t`, `x`, `y`
#' (mean position per window).
#' @param x x positions of one animal.
#' @param t time in seconds, 0 at lights on.
#' @param time_window_length window size in seconds.
#' @param context_bins windows on each side used to judge stillness (3, i.e. 30 s).
#' @param max_shift largest positional change still counted as still, in the
#' units of `x`/`y`. `NULL` uses one pixel ([default_still_shift]).
#' @param day_length length of the day in hours.
#' @param lights_off hour of lights off.
#' @seealso [max_velocity_detector], [motion_qc], [sleep_rules]
NULL

#' @describeIn motion_calibration one pixel in the units of `x`: positions
#' normalised to the ROI width (as scopr loads them) or in pixels.
#' @export
pixel_size <- function(x){
  x <- x[is.finite(x)]
  # positions normalised to ROI width lie in [0, 1]; 2e-3 of a ~500 px ROI is 1 px
  if(length(x) > 0 && max(x) <= 1.5)
    return(2e-3)
  1
}

#' @describeIn motion_calibration positional tolerance for a still window, 1 pixel in
#' the units of `x`. Validated against tracker-independent pixel motion on two
#' recordings: at 1 px the calibrated threshold matched the one chosen from the
#' pixels; at 0.5 px it skipped still flies whose tracked position jitters and came
#' out too low.
#' @export
default_still_shift <- function(x){
  STILL_SHIFT_PIXELS <- 1
  STILL_SHIFT_PIXELS * pixel_size(x)
}

#' @describeIn motion_calibration `"light"` / `"dark"` label for each time point.
#' @export
light_phase <- function(t, day_length = 24, lights_off = 12){
  ifelse((t %% (day_length * 3600)) < lights_off * 3600, "light", "dark")
}

#' @describeIn motion_calibration tracking spikes: the centroid jumps more than
#' `min_jump` away for up to `max_frames` frames and lands back within
#' `tolerance` of where it was, within `max_span` seconds -- the tracker briefly
#' locking onto something other than the animal. Returns a list of two logical
#' vectors over frames: `velocity` (displaced frames and the return frame) and
#' `position` (displaced frames only).
#' @param y y positions.
#' @param max_frames longest displaced run treated as a spike.
#' @param min_jump smallest displacement treated as a jump; `NULL` is 3 pixels.
#' @param tolerance how close the return must be; `NULL` is 1 pixel.
#' @param max_span longest time from the last good frame to the return, in seconds.
#' @export
find_spikes <- function(t, x, y, max_frames = 2, min_jump = NULL,
                        tolerance = NULL, max_span = 2){
  px <- pixel_size(x)
  if(is.null(min_jump)) min_jump <- 3 * px
  if(is.null(tolerance)) tolerance <- px
  n <- length(t)
  velocity <- rep(FALSE, n)
  position <- rep(FALSE, n)
  for(m in seq_len(max_frames)){
    if(n < m + 2)
      break
    # anchor a = last good frame, displaced a+1 .. a+m, return a+m+1
    a <- seq_len(n - m - 1)
    displaced <- rep(TRUE, length(a))
    for(j in seq_len(m))
      displaced <- displaced & sqrt((x[a + j] - x[a])^2 + (y[a + j] - y[a])^2) > min_jump
    back <- sqrt((x[a + m + 1] - x[a])^2 + (y[a + m + 1] - y[a])^2) <= tolerance
    quick <- (t[a + m + 1] - t[a]) <= max_span
    hit <- a[which(displaced & back & quick)]
    for(j in seq_len(m)){
      velocity[hit + j] <- TRUE
      position[hit + j] <- TRUE
    }
    velocity[hit + m + 1] <- TRUE
  }
  list(velocity = velocity, position = position)
}

#' Quantile of each row of a matrix, ignoring NA
#'
#' Same definition as `stats::quantile(type = 7)` (and numpy's default), but
#' vectorised over rows: `apply` is ~100x slower on the (n, 7) windows used by
#' [find_still_bins].
#' @noRd
row_quantile <- function(m, q){
  nr <- nrow(m)
  nc <- ncol(m)
  # sort within rows in one call; NA go last within each row
  o <- order(row(m), m, na.last = TRUE)
  ordered <- matrix(m[o], nrow = nr, ncol = nc, byrow = TRUE)
  count <- rowSums(is.finite(m))
  pos <- pmax(count - 1, 0) * q
  lower <- floor(pos)
  upper <- ceiling(pos)
  low_val <- ordered[cbind(seq_len(nr), lower + 1)]
  high_val <- ordered[cbind(seq_len(nr), upper + 1)]
  out <- low_val + (high_val - low_val) * (pos - lower)
  out[count == 0] <- NA_real_
  out
}

#' @describeIn motion_calibration logical vector, `TRUE` for windows in which the
#' animal did not change position. The median positions before and after the
#' window must agree, and the positions across the whole window must agree once
#' the single highest and lowest window are set aside, so that a glitch in one
#' window is forgiven but pacing back and forth is not.
#' @export
find_still_bins <- function(d_small,
                            time_window_length = 10,
                            context_bins = 3,
                            max_shift = NULL){
  if(is.null(max_shift))
    max_shift <- default_still_shift(d_small$x)

  grid <- seq(from = min(d_small$t), to = max(d_small$t), by = time_window_length)
  # regular grid so windows count time windows, not rows, across tracking gaps
  idx <- match(d_small$t, grid)
  width <- 2 * context_bins + 1
  still <- rep(FALSE, length(grid))
  if(length(grid) >= width){
    n_win <- length(grid) - width + 1
    min_side <- max(1, context_bins - 1)
    ok <- rep(TRUE, n_win)
    shift_sq <- rep(0, n_win)
    before_cols <- seq_len(context_bins)
    after_cols <- (context_bins + 2):width
    for(axis in c("x", "y")){
      v <- rep(NA_real_, length(grid))
      v[idx] <- d_small[[axis]]
      win <- stats::embed(v, width)[, width:1, drop = FALSE]
      n_before <- rowSums(is.finite(win[, before_cols, drop = FALSE]))
      n_after <- rowSums(is.finite(win[, after_cols, drop = FALSE]))
      ok <- ok & n_before >= min_side & n_after >= min_side
      shift_sq <- shift_sq + (row_quantile(win[, after_cols, drop = FALSE], 0.5) -
                                row_quantile(win[, before_cols, drop = FALSE], 0.5))^2
      # trimming one window at each end forgives a single glitched window,
      # while back-and-forth movement spreads over several windows and is kept
      spread <- row_quantile(win, 5/6) - row_quantile(win, 1/6)
      ok <- ok & !is.na(spread) & spread <= max_shift
    }
    ok <- ok & !is.na(shift_sq) & sqrt(shift_sq) <= max_shift
    still[(context_bins + 1):(context_bins + n_win)] <- ok
  }
  still[idx]
}

#' Check how well the fixed movement threshold fits the tracking noise
#'
#' Run on raw ethoscope data before scoring sleep. For each animal and light
#' phase it reports how often a still animal crosses the fixed threshold
#' (`fp_rate_fixed`) and what fraction of genuine rests of `min_time_immobile`
#' would survive that (`rest_survival_fixed`, `(1 - fp)^windows`).
#' `fp_rate_fixed` above about 0.01 means classic sleep estimates are unreliable;
#' score sleep with `sleep_annotation(rule = "k")` instead. A large
#' `untracked_fraction` matters under either rule: with the default
#' `untracked = "immobile"`, windows where the animal was lost count as still.
#'
#' @param data [behavr::behavr] table (or keyed [data.table::data.table]) of raw
#' ethoscope data, or the data of a single animal.
#' @param velocity_threshold the fixed threshold being checked.
#' @param min_time_immobile sleep criterion in seconds.
#' @inheritParams motion_calibration
#' @inheritParams motion_detectors
#' @return a [data.table::data.table] with one row per animal and phase:
#' `phase`, `n_bins`, `untracked_fraction`, `spike_fraction` (frames that are
#' tracking spikes, see [find_spikes]), `n_still`, `fp_rate_fixed`,
#' `rest_survival_fixed`.
#' @seealso [motion_calibration], [sleep_rules]
#' @export
motion_qc <- function(data,
                      time_window_length = 10,
                      velocity_correction_coef = 3e-3,
                      velocity_threshold = 1,
                      min_time_immobile = 300,
                      day_length = 24,
                      lights_off = 12){
  .SD = x = y = t = NULL
  wrapped <- function(d){
    if(nrow(d) < 100)
      return(NULL)
    d_small <- max_velocity_detector(d, time_window_length,
                                     velocity_correction_coef = velocity_correction_coef,
                                     masking_duration = 0)
    if(!"y" %in% colnames(d))
      d <- data.table::copy(d)[, y := 0]
    dd <- prepare_data_for_motion_detector(d, c("t", "x", "y"), time_window_length)
    pos <- dd[, list(x = mean(x), y = mean(y)), by = "t_round"]
    binned <- data.table::data.table(t = d_small$t, x = pos$x, y = pos$y,
                                     max_velocity = d_small$max_velocity)
    frames <- d[order(t)]
    spikes <- find_spikes(frames$t, frames$x, frames$y)$velocity
    frame_phase <- light_phase(frames$t, day_length, lights_off)
    still <- find_still_bins(binned, time_window_length)
    phase <- light_phase(binned$t, day_length, lights_off)
    grid <- seq(from = min(binned$t), to = max(binned$t), by = time_window_length)
    grid_phase <- light_phase(grid, day_length, lights_off)
    out <- lapply(intersect(c("light", "dark"), unique(phase)), function(name){
      in_phase <- phase == name
      mask <- still & in_phase & is.finite(binned$max_velocity)
      fp <- if(any(mask)) mean(binned$max_velocity[mask] > velocity_threshold) else NA_real_
      data.table::data.table(phase = name,
                             n_bins = sum(in_phase),
                             untracked_fraction = 1 - sum(in_phase) / sum(grid_phase == name),
                             spike_fraction = mean(spikes[frame_phase == name]),
                             n_still = sum(mask),
                             fp_rate_fixed = fp,
                             rest_survival_fixed = (1 - fp)^(min_time_immobile / time_window_length))
    })
    data.table::rbindlist(out)
  }
  if(is.null(data.table::key(data)))
    return(wrapped(data))
  data[, wrapped(.SD), by = eval(data.table::key(data))]
}
