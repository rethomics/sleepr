#' Sleep from walking and sustained movement (the k-rule)
#'
#' The classic rule calls a 10-s window "moving" when any frame's velocity
#' exceeds a threshold, so a single noisy frame breaks a sleep bout. The k-rule
#' asks two questions instead: did the animal walk (its median position moved
#' more than 10 px from the previous window), and did it move in place (at least
#' `k` sustained movement events started in a centred 60-s window)?
#'
#' A movement event is a run of consecutive frames above the classic velocity
#' test. Events in which the position never leaves the pixel it started from
#' ("subpixel"), and those that jump out and land back within a pixel in one or
#' two frames ("flicker"), are tracking noise and are ignored; only the remaining
#' "sustained" events count. Sleep is then 5 minutes or more of tracked windows
#' with neither walking nor such movement. A window without frames is never sleep.
#'
#' The rule is tentative and opt-in, through `sleep_annotation(rule = "k")`. The
#' algorithm and its defaults are identical to ethoscopy's `sleep_rules` module,
#' so both packages score the same data alike.
#'
#' @name sleep_rules
#' @param x,y positions of one animal, one per frame, in recording order.
#' @param moving logical, `TRUE` where the frame passes the movement test.
#' @param tolerance how close a flicker must return, in the units of `x`/`y`.
#' @param max_flicker_frames longest run of moving frames that can be a flicker.
#' @param t frame times in seconds, in recording order.
#' @param xy_dist_log10x1000 the tracker's per-frame displacement.
#' @param k sustained events in the centred 60-s window that make a window awake.
#' @param pixel one pixel in the units of `x`/`y`.
#' @param velocity_correction_coef as in [max_velocity_detector]; a frame moves
#' when its corrected velocity exceeds 1.
#' @param min_sleep_bins shortest sleep bout, in 10-s windows.
#' @seealso [sleep_annotation]
NULL

K_RULE_BIN_SECONDS <- 10
WALK_SHIFT_PIXELS <- 10
RETURN_TOLERANCE_PIXELS <- 1
MAX_FLICKER_FRAMES <- 2
# events are counted in windows [i - 3, i + 3)
EVENT_WINDOW_BINS <- 6

#' @describeIn sleep_rules split runs of moving frames into events. The frame
#' before a run is its anchor, so a run on the first frame is not an event. An
#' event is subpixel if every frame is exactly at the anchor's position, else a
#' flicker if it has at most `max_flicker_frames` frames and the next frame is
#' back within `tolerance` of the anchor, else sustained. Returns a list with
#' `start` (row of the first frame) and `class` (0 subpixel, 1 flicker, 2 sustained).
#' @export
classify_events <- function(x, y, moving,
                            tolerance = RETURN_TOLERANCE_PIXELS,
                            max_flicker_frames = MAX_FLICKER_FRAMES){
  x <- as.numeric(x)
  y <- as.numeric(y)
  moving <- as.logical(moving)
  moving[is.na(moving)] <- FALSE
  runs <- rle(moving)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1
  keep <- runs$values & starts > 1
  starts <- starts[keep]
  ends <- ends[keep]
  if(length(starts) == 0)
    return(list(start = integer(0), class = integer(0)))
  anchor <- starts - 1
  lengths <- ends - starts + 1
  run_of_frame <- rep(seq_along(starts), lengths)
  frame <- sequence(lengths, from = starts)
  away <- x[frame] != x[anchor[run_of_frame]] | y[frame] != y[anchor[run_of_frame]]
  # NaN positions differ from everything, as in numpy
  away[is.na(away)] <- TRUE
  away_count <- c(0, cumsum(away))
  subpixel <- away_count[cumsum(lengths) + 1] - away_count[cumsum(lengths) - lengths + 1] == 0
  has_next <- ends < length(x)
  after <- pmin(ends + 1, length(x))
  back <- sqrt((x[after] - x[anchor])^2 + (y[after] - y[anchor])^2) <= tolerance
  back[is.na(back)] <- FALSE
  flicker <- lengths <= max_flicker_frames & has_next & back
  list(start = starts,
       class = ifelse(subpixel, 0L, ifelse(flicker, 1L, 2L)))
}

#' @describeIn sleep_rules score one animal's frames window by window. Returns a
#' [data.table::data.table] with one row per 10-s window from the first to the
#' last window with frames: `t` (window start, s), `has_data`, `walking`,
#' `sustained` (events starting in the window), `micro_awake` and `asleep`.
#' @export
k_rule_bins <- function(t, x, y, xy_dist_log10x1000,
                        k = 3,
                        pixel = 1,
                        velocity_correction_coef = 3e-3,
                        min_sleep_bins = 30){
  bin = n = .N = . = NULL
  if(length(t) == 0)
    return(data.table::data.table(t = numeric(0), has_data = logical(0),
                                  walking = logical(0), sustained = integer(0),
                                  micro_awake = logical(0), asleep = logical(0)))
  x <- as.numeric(x)
  y <- as.numeric(y)
  # floor() rather than %/%, which can be one less for non-integer divisions
  frame_bin <- floor(as.numeric(t) / K_RULE_BIN_SECONDS)
  frames <- data.table::data.table(bin = frame_bin, x = x, y = y)
  medians <- frames[, .(x = stats::median(x, na.rm = TRUE),
                        y = stats::median(y, na.rm = TRUE),
                        n = .N), by = bin]
  first <- min(frame_bin)
  grid <- medians[data.table::data.table(bin = seq(first, max(frame_bin))), on = "bin"]
  has_data <- !is.na(grid$n)
  step <- sqrt(c(NA, diff(grid$x))^2 + c(NA, diff(grid$y))^2)
  # an unknown step (first window, or a neighbouring window without frames) counts as walking
  still <- has_data & !is.na(step) & step <= WALK_SHIFT_PIXELS * pixel

  moving <- 10 ^ (as.numeric(xy_dist_log10x1000) / 1000) / velocity_correction_coef > 1
  events <- classify_events(x, y, moving, tolerance = RETURN_TOLERANCE_PIXELS * pixel)
  n_bins <- nrow(grid)
  event_bin <- frame_bin[events$start[events$class == 2L]] - first + 1
  sustained <- tabulate(event_bin, nbins = n_bins)
  cumulative <- c(0, cumsum(sustained))
  half <- EVENT_WINDOW_BINS %/% 2
  i <- seq_len(n_bins) - 1
  lo <- pmin(pmax(i - half, 0), n_bins)
  hi <- pmin(pmax(i + EVENT_WINDOW_BINS - half, 0), n_bins)
  micro_awake <- cumulative[hi + 1] - cumulative[lo + 1] >= k

  runs <- rle(still & !micro_awake)
  asleep <- rep(runs$values & runs$lengths >= min_sleep_bins, runs$lengths)
  data.table::data.table(t = grid$bin * K_RULE_BIN_SECONDS,
                         has_data = has_data,
                         walking = has_data & !still,
                         sustained = sustained,
                         micro_awake = micro_awake,
                         asleep = asleep)
}

#' frames the tracker observed: `is_inferred` numerically 0 (it can be text)
#' @noRd
observed_frames <- function(is_inferred){
  if(is.logical(is_inferred) || is.numeric(is_inferred))
    value <- as.numeric(is_inferred)
  else
    value <- suppressWarnings(as.numeric(as.character(is_inferred)))
  !is.na(value) & value == 0
}

#' argument checks for sleep_annotation(rule = "k")
#' @noRd
check_k_rule_args <- function(time_window_length, k, pixel, dots){
  if("velocity_threshold" %in% names(dots))
    stop('velocity_threshold applies only to rule = "classic"')
  if(time_window_length != K_RULE_BIN_SECONDS)
    stop(sprintf('rule = "k" is defined for %d-s windows', K_RULE_BIN_SECONDS))
  if(!is.numeric(k) || length(k) != 1 || is.na(k) || k < 1 || k != round(k))
    stop("k must be a positive integer")
  if(!is.null(pixel) && !(is.numeric(pixel) && length(pixel) == 1 && !is.na(pixel) && pixel > 0))
    stop("pixel must be positive")
}

#' the body of sleep_annotation(rule = "k") for one animal
#' @noRd
k_rule_annotation <- function(d, time_window_length, min_time_immobile,
                              motion_detector_FUN, k, pixel, columns_to_keep, ...){
  moving = is_interpolated = walking = sustained = micro_awake = asleep = t = NULL
  if("is_inferred" %in% colnames(d))
    d <- d[observed_frames(d$is_inferred)]
  if(nrow(d) < 100)
    return(NULL)
  dots <- list(...)
  coef <- if(is.null(dots$velocity_correction_coef)) 3e-3 else dots$velocity_correction_coef
  kb <- k_rule_bins(d$t, d$x, d$y, d$xy_dist_log10x1000,
                    k = k,
                    pixel = if(is.null(pixel)) pixel_size(d$x) else pixel,
                    velocity_correction_coef = coef,
                    min_sleep_bins = min_time_immobile / time_window_length)
  time_map <- data.table::data.table(t = kb$t, key = "t")
  d_small <- motion_detector_FUN(d, time_window_length, ...)
  # the rule is scored even where the classic detector found too few frames
  if(is.null(d_small) || nrow(d_small) < 1){
    d_small <- data.table::copy(time_map)[, moving := FALSE]
  } else {
    missing_val <- time_map[!d_small]
    d_small <- d_small[time_map, roll = TRUE]
    d_small[t %in% missing_val$t, moving := FALSE]
  }
  d_small[, `:=`(is_interpolated = !kb$has_data,
                 walking = kb$walking,
                 sustained = kb$sustained,
                 micro_awake = kb$micro_awake,
                 asleep = kb$asleep)]
  # keep every window of the rule's grid, even where frame-level columns are missing
  d_small <- stats::na.omit(d[d_small, on = c("t"), roll = TRUE], cols = c("t", "asleep"))
  d_small[, intersect(columns_to_keep, colnames(d_small)), with = FALSE]
}
