context("motion_calibration")

# Synthetic fly with known truth: rests at a fixed position and walks for the
# first `active_minutes` of every hour. Positions are normalised to the ROI width,
# as scopr loads them. A 2-hour "day" (1 h light, 1 h dark) keeps data small.
roi_width <- 508
make_track <- function(hours = 4, fps = 2, noise_light = 0.9, noise_dark = 0.9,
                       active_minutes = 10, always_active = FALSE, gap = NULL,
                       spike_rate = 0, seed = 0){
  set.seed(seed)
  t <- seq(0, hours * 3600 - 1 / fps, by = 1 / fps)
  active <- always_active | ((t %% 3600) < active_minutes * 60)
  x_px <- 200 + ifelse(active, 60 * sin(2 * pi * t / 20), 0)
  x_px <- round(x_px)
  # tracker jumps 40 px for one frame and lands back on the same pixel
  x_px <- x_px + ifelse(!active & runif(length(t)) < spike_rate, 40, 0)
  noise <- ifelse(light_phase(t, 2, 1) == "light", noise_light, noise_dark)
  # a still fly's distance is lognormal jitter around the noise level
  jitter <- noise * exp(rnorm(length(t), 0, 0.08)) * 3e-3
  step <- abs(c(0, diff(x_px))) / roi_width
  d <- data.table::data.table(t = t, x = x_px / roi_width, y = 30 / roi_width,
                              xy_dist_log10x1000 = round(1000 * log10(step + jitter)),
                              has_interacted = 0L)
  if(!is.null(gap))
    d <- d[t < gap[1] | t >= gap[2]]
  d
}

rest_sleep <- function(res, active_minutes = 10, settle_minutes = 6)
  res[(t %% 3600) >= (active_minutes + settle_minutes) * 60, mean(asleep)]

test_that("invalid arguments are rejected", {
  d <- make_track(hours = 1)
  expect_error(sleep_annotation(d, rule = "classic", untracked = "skip"))
})

test_that("untracked = 'break' ends sleep at untracked windows", {
  gap <- c(3600 + 30 * 60, 3600 + 40 * 60)
  d <- make_track(noise_light = 0.3, noise_dark = 0.3, gap = gap)
  immobile <- sleep_annotation(d, rule = "classic", untracked = "immobile")
  broken <- sleep_annotation(d, rule = "classic", untracked = "break")
  in_gap <- immobile$t >= gap[1] & immobile$t < gap[2]
  expect_true(all(immobile$is_interpolated[in_gap]))
  expect_true(all(immobile$asleep[in_gap]))
  expect_false(any(broken$asleep[broken$t >= gap[1] & broken$t < gap[2]]))
  expect_identical(immobile$moving, broken$moving)
})

test_that("still bins forgive a glitch but not a shift or pacing", {
  binned <- function(x) data.table::data.table(t = (seq_along(x) - 1) * 10, x = x, y = 30)
  glitch <- rep(200, 40); glitch[21] <- 203
  expect_true(all(find_still_bins(binned(glitch))[6:35]))
  shift <- find_still_bins(binned(c(rep(200, 20), rep(205, 20))))
  expect_false(shift[20] || shift[21])
  expect_true(shift[6] && shift[36])
  # sub-pixel jitter of a still fly is still at the default 1 px, not at 0.5 px
  jitter <- 200 + rep(c(0, 0.4, 0.8, 0.3, 0.6), 8)
  expect_true(all(find_still_bins(binned(jitter))[6:35]))
  expect_lt(mean(find_still_bins(binned(jitter), max_shift = 0.5)[6:35]), 0.8)
  pacing <- find_still_bins(binned(rep(c(162, 238), 20)))
  expect_false(any(pacing))
  edges <- find_still_bins(binned(rep(200, 20)))
  expect_false(edges[1] || edges[20])
})

test_that("default still shift follows position units", {
  expect_equal(default_still_shift(c(12, 480)), 1)
  expect_equal(default_still_shift(c(0.02, 0.95)), 2e-3)
})

test_that("motion_qc flags a noisy recording", {
  noisy <- make_track()[, id := "noisy"]
  clean <- make_track(noise_light = 0.3, noise_dark = 0.3, seed = 1)[, id := "clean"]
  d <- rbind(noisy, clean)
  data.table::setkeyv(d, "id")
  qc <- motion_qc(d, day_length = 2, lights_off = 1)
  expect_equal(nrow(qc), 4)
  expect_true(all(qc[id == "noisy", fp_rate_fixed] > 0.5))
  expect_true(all(qc[id == "noisy", rest_survival_fixed] < 0.01))
  expect_true(all(qc[id == "clean", fp_rate_fixed] == 0))
  expect_false("auto_threshold" %in% names(qc))
})

test_that("still bins match ethoscopy on the same binned data", {
  # parity_binned.csv: binned data and ethoscopy's still bins
  ref <- data.table::fread(test_path("parity_binned.csv"))
  expect_identical(find_still_bins(ref), as.logical(ref$still))
})

test_that("one and two frame spikes are found, real moves are not", {
  tt <- (0:29) * 0.25
  x <- rep(200, 30); x[11] <- 240; x[21:22] <- 230
  sp <- find_spikes(tt, x, rep(30, 30))
  expect_equal(which(sp$position), c(11, 21, 22))
  expect_equal(which(sp$velocity), c(11, 12, 21, 22, 23))
  step <- c(rep(200, 15), rep(240, 15))
  jitter <- rep(c(200, 201), 15)
  for(v in list(step, jitter))
    expect_false(any(find_spikes(tt, v, rep(30, 30))$velocity))
  far <- find_spikes(c(0, .25, .5, 10, 10.25), c(200, 200, 240, 200, 200), rep(30, 5))
  expect_false(any(far$velocity))
})

test_that("motion_qc reports spikes", {
  spiky <- make_track(spike_rate = 0.004)
  qc <- motion_qc(spiky, day_length = 2, lights_off = 1)
  expect_true(all(qc$spike_fraction > 0.002))
})

test_that("spike detection matches ethoscopy on the same frames", {
  ref <- data.table::fread(test_path("parity_frames.csv"))
  sp <- find_spikes(ref$t, ref$x, ref$y)
  expect_identical(sp$velocity, as.logical(ref$spike_velocity))
  expect_identical(sp$position, as.logical(ref$spike_position))
})
