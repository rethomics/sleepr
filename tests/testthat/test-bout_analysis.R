context("bout_analysis")

test_that("bout_analysis works resturn expected results", {
  library(behavr)
  dt <- toy_activity_data(duration=days(2))
  bout_analysis(asleep, dt)
})

