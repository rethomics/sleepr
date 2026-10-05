context("sleep_rules")

STILL <- -3000  # xy_dist_log10x1000 for a corrected velocity of ~0.33
MOVE <- -2000   # ~3.3

classes <- function(x, moving, y = rep(0, length(x))){
  ev <- classify_events(x, y, as.logical(moving))
  list(start = ev$start, class = ev$class)
}

# one x position per 10-s window, 10 frames per window at 1 fps
bins <- function(x, d = NULL, shift = NULL, k = 3, ...){
  x <- rep(x, each = 10)
  xy <- rep(STILL, length(x))
  if(!is.null(d)) xy[d] <- MOVE
  if(!is.null(shift)) x[shift] <- 130
  t <- seq_along(x) - 1
  keep <- !is.na(x)
  k_rule_bins(t[keep], x[keep], rep(0, sum(keep)), xy[keep], k = k, ...)
}

# one sustained event per window: a moving frame that leaves its pixel and does
# not come back on the next frame (frames are 1-based, as in R)
sustained_at <- function(windows){
  frames <- windows * 10 + 4
  list(d = frames, shift = c(frames, frames + 1))
}

test_that("events are subpixel, flicker or sustained", {
  ev <- classes(c(5, 5, 5, 9, 5, 5, 8, 9, 5, 5, 7, 8, 9, 5),
                c(0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0))
  expect_equal(ev$start, c(2, 4, 7, 11))
  expect_equal(ev$class, c(0L, 1L, 1L, 2L))
})

test_that("the return tolerance is inclusive", {
  expect_equal(classes(c(0, 4, 1), c(0, 1, 0))$class, 1L)
  expect_equal(classes(c(0, 4, 1.01), c(0, 1, 0))$class, 2L)
  expect_equal(classes(c(0, 4, 0.6), c(0, 1, 0), y = c(0, 4, 0.8))$class, 1L)
})

test_that("runs on the first and last frame", {
  expect_equal(classes(c(5, 9, 5, 9, 5), c(1, 1, 0, 1, 0)), list(start = 4, class = 1L))
  expect_equal(classes(c(5, 5, 9), c(0, 0, 1)), list(start = 3, class = 2L))
  expect_equal(classes(c(5, 5, 5), c(0, 0, 1)), list(start = 3, class = 0L))
  expect_length(classify_events(1:5, rep(0, 5), rep(FALSE, 5))$start, 0)
})

test_that("walking threshold, bouts and windows without frames", {
  expect_equal(bins(c(100, 110, 120.5, 120.5))$walking, c(TRUE, FALSE, TRUE, FALSE))
  expect_false(any(bins(rep(100, 30))$asleep))
  expect_equal(sum(bins(rep(100, 31))$asleep), 30)
  out <- bins(c(rep(100, 35), NA, rep(100, 35)))
  expect_false(out$has_data[36])
  expect_true(out$walking[37])
  expect_true(all(out$asleep[c(2:35, 38:71)]))
  expect_false(any(out$asleep[c(1, 36, 37)]))
})

test_that("sustained events, the 60-s window and k", {
  ev <- sustained_at(20:22)
  out <- bins(rep(100, 60), d = ev$d, shift = ev$shift, k = 3)
  expect_equal(sum(out$sustained), 3)
  expect_equal(which(out$micro_awake) - 1, 20:23)
  two <- bins(rep(100, 60), d = ev$d, shift = ev$shift, k = 2)
  expect_equal(which(two$micro_awake) - 1, 19:24)
  ends <- sustained_at(c(0, 39))
  clipped <- bins(rep(100, 40), d = ends$d, shift = ends$shift, k = 1)
  expect_equal(which(clipped$micro_awake) - 1, c(0:3, 37:39))
})

test_that("pixel scales the limits", {
  x <- c(100, 109, 125, 125, 125)
  expect_equal(bins(x)$walking, c(TRUE, FALSE, TRUE, FALSE, FALSE))
  expect_equal(bins(x / 551, pixel = 1 / 551)$walking, bins(x)$walking)
})

raw_track <- function(n_bins = 120, inferred = NULL){
  n <- n_bins * 10
  d <- data.table::data.table(t = seq_len(n) - 1, x = 100, y = 20,
                              xy_dist_log10x1000 = STILL, has_interacted = 0L)
  if(!is.null(inferred)) d[, is_inferred := inferred]
  d
}

test_that("sleep_annotation(rule = 'k') drops inferred frames and keeps the grid", {
  inferred <- ifelse(seq_len(1600) > 800 & seq_len(1600) <= 820, "1", "0")
  out <- sleep_annotation(raw_track(160, inferred), rule = "k", masking_duration = 0)
  expect_true(all(c("walking", "sustained", "micro_awake", "is_interpolated") %in% names(out)))
  expect_equal(out$t, (0:159) * 10)
  expect_true(all(out$is_interpolated[81:82]))
  expect_true(out$walking[83])
  expect_true(all(out$asleep[c(2:80, 84:160)]))
  expect_false(any(out$asleep[c(1, 81:83)]))
})

test_that("sleep_annotation(rule = 'k') rejects what the rule does not take", {
  d <- raw_track()
  expect_error(sleep_annotation(d, rule = "k", velocity_threshold = 1))
  expect_error(sleep_annotation(d, rule = "k", time_window_length = 20))
  expect_error(sleep_annotation(d, rule = "k", k = 0))
  expect_error(sleep_annotation(d, rule = "k", k = TRUE))
  expect_error(sleep_annotation(d, rule = "k", k = 2.5))
  expect_error(sleep_annotation(d, rule = "k", pixel = 0))
  expect_error(sleep_annotation(d, rule = "nope"))
  expect_null(sleep_annotation(raw_track(5), rule = "k"))
})

test_that("the k-rule asks scopr for y positions", {
  needed <- attr(sleep_annotation, "needed_columns")
  expect_false("y" %in% needed())
  expect_true(all(c("x", "y", "xy_dist_log10x1000") %in% needed(rule = "k")))
})

test_that("k-rule matches ethoscopy and the reference bin by bin", {
  # k_rule_*.csv: real recordings in pixels, exported from ethoscopy with the
  # reference implementation's per-window results (scripts/validate_k_rule.py)
  frames <- data.table::fread(test_path("k_rule_frames.csv"))
  expected <- data.table::fread(test_path("k_rule_bins.csv"))
  expect_true(length(unique(frames$fly)) >= 3)
  for(fly_id in unique(frames$fly)){
    f <- frames[frames$fly == fly_id]
    f <- f[observed_frames(f$is_inferred)]
    ref <- expected[expected$fly == fly_id]
    for(k in c(3, 2)){
      out <- k_rule_bins(f$t / 1000, f$x, f$y, f$xy_dist_log10x1000, k = k, pixel = 1)
      expect_equal(out$t, ref$t, info = fly_id)
      expect_identical(out$has_data, ref$has_data, info = fly_id)
      expect_identical(out$walking, ref$walking, info = fly_id)
      expect_identical(out$asleep, ref[[paste0("asleep_k", k)]], info = paste(fly_id, k))
    }
  }
})

test_that("sleep_annotation(rule = 'k') reproduces the reference on a real recording", {
  frames <- data.table::fread(test_path("k_rule_frames.csv"))
  expected <- data.table::fread(test_path("k_rule_bins.csv"))
  f <- frames[frames$fly == "abg_2026"][, `:=`(t = t / 1000, has_interacted = 0L, fly = NULL)]
  out <- sleep_annotation(f, rule = "k", pixel = 1)
  ref <- expected[expected$fly == "abg_2026"]
  expect_equal(out$t, ref$t)
  expect_identical(out$asleep, ref$asleep_k3)
  expect_identical(!out$is_interpolated, ref$has_data)
})
