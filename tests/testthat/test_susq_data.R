# tests/testthat/test_susq_data.R
source("../../R/susq_data.R")

test_that("prepare_susq_data returns 164 determinations from 167 rows", {
  d <- prepare_susq_data()
  expect_equal(nrow(d$dates), 164L)
  expect_equal(length(unique(d$dates$site)), 34L)
})

test_that("the failed Bates_66 replicate is pooled by the Ward and Wilson multiplier", {
  d <- prepare_susq_data()
  # UGAMS-53046 (546 +/- 20) and UGAMS-53046r (637 +/- 25): T = 8.08 on 1 df.
  # Pooled mean 581.5, pooled error 15.61, multiplier sqrt(8.08) = 2.843.
  row <- d$dates[d$dates$lab_no == "UGAMS-53046", ]
  expect_equal(nrow(row), 1L)
  expect_equal(row$c14_age, 581.5, tolerance = 0.1)
  expect_equal(row$c14_error, 44.4, tolerance = 0.5)
  expect_equal(row$preprocess, "pooled_replicate")
  expect_false(any(d$dates$lab_no == "UGAMS-53046r"))
})

test_that("wiggle-matched pairs become single determinations carrying their ring offset", {
  d <- prepare_susq_data()
  wm <- d$dates[d$dates$preprocess == "wiggle_match", ]
  expect_equal(nrow(wm), 2L)
  b46 <- wm[wm$lab_no == "UGAMS-59365", ]
  expect_equal(b46$c14_age, 870)
  expect_equal(b46$pair_lab_no, "UGAMS-59366")
  expect_equal(b46$pair_c14_age, 820)
  expect_equal(b46$pair_offset, 17)   # RY3 to RY20 midpoints
  b89 <- wm[wm$lab_no == "UGAMS-19367", ]
  expect_equal(b89$pair_offset, 25)   # RY3 to RY28 midpoints
})

test_that("material classes match the audited lookup counts", {
  d <- prepare_susq_data()
  # 78 short_lived, 68 wood, 21 indeterminate across 167 rows; after
  # preprocessing, 3 rows are absorbed: one replicate and two pair members,
  # all of them wood or short_lived.
  tab <- table(d$dates$material_class)
  expect_equal(as.integer(tab[["short_lived"]]), 77L)  # one maize replicate absorbed
  expect_equal(as.integer(tab[["wood"]]), 66L)         # two ring-pair members absorbed
  expect_equal(as.integer(tab[["indeterminate"]]), 21L)
  expect_true(all(d$dates$class_id == match(d$dates$material_class,
                  c("short_lived", "wood", "indeterminate"))))
})

test_that("phase counts follow the 6 / 3 rule and total 74", {
  d <- prepare_susq_data()
  pc <- site_phase_counts(d$dates)
  expect_equal(sum(pc$J), 74L)
  expect_equal(sum(pc$J == 4L), 9L)
  expect_equal(sum(pc$J == 2L), 13L)
  expect_equal(sum(pc$J == 1L), 12L)
  expect_true(all(pc$J[pc$n >= 6] == 4L))
  expect_true(all(pc$J[pc$n <= 2] == 1L))
})

test_that("quality flags are reported and not silently corrected", {
  d <- prepare_susq_data()
  expect_true(any(grepl("AA-41933", d$quality$detail)))
  expect_true(any(grepl("none provided", d$quality$detail)))
  expect_true(any(grepl("Beta-7007", d$quality$detail)))
})

test_that("pooling a replicate is itself reported, with the test that justified it", {
  # Pooling alters a determination. The spec requires data problems to be
  # reported, so the inconsistency that forced the error inflation must appear
  # in the quality report rather than only inside the pooling function.
  d <- prepare_susq_data()
  row <- d$quality[d$quality$kind == "replicate_inconsistent", ]
  expect_equal(nrow(row), 1L)
  expect_match(row$detail, "UGAMS-53046")
  expect_match(row$detail, "T = 8.08")
  expect_match(row$detail, "not consistent")
})
