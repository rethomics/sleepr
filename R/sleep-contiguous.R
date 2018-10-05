#' @param moving logical vector of acivity
#' @param fs sampling frequency (Hz)
#' @param min_valid_time the minimal amount immobile time that counts as sleep (s)
#' @return a boolean vector of the same length as `moving`.
#' `TRUE` values where and only where `!moving` for a run length of at least 5min
#' @seealso [rle] on which this function is based
#' @noRd
sleep_contiguous <- function(moving, fs, min_valid_time = 5*60){
  min_len <- fs * min_valid_time
  r_sleep <- rle(!moving)
  valid_runs <-  r_sleep$length >= min_len
  r_sleep$values <- valid_runs & r_sleep$value
  inverse.rle(r_sleep)
}


