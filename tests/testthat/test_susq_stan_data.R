source("../../R/susq_data.R")
source("../../R/susq_calibration.R")
source("../../R/susq_stan_data.R")

d <- prepare_susq_data()
cc <- intcal20_on_grid()
v <- build_variants(d$dates, cc)

test_that("short-lived determinations get 2 variants and others get 7", {
  n_by_det <- table(v$index$det_id)
  cls <- d$dates$material_class[match(as.integer(names(n_by_det)), d$dates$det_id)]
  expect_true(all(n_by_det[cls == "short_lived"] == 2L))
  expect_true(all(n_by_det[cls != "short_lived"] == 7L))
  expect_equal(nrow(v$index), 77L * 2L + (66L + 21L) * 7L)
})

test_that("every determination has exactly one outlier variant and one kappa=0 variant", {
  by_det <- split(v$index, v$index$det_id)
  expect_true(all(vapply(by_det, function(x) sum(x$is_outlier), 0) == 1))
  expect_true(all(vapply(by_det, function(x) sum(x$kappa_idx == 1L), 0) == 1))
})

test_that("variants of one determination share a normaliser, so relative scale survives", {
  # The kappa = 0 variant integrates to 1 by construction. The outlier variant
  # must NOT also integrate to 1, or the outlier evidence has been destroyed.
  i <- v$index$variant_id[v$index$det_id == 1 & v$index$kappa_idx == 1L]
  o <- v$index$variant_id[v$index$det_id == 1 & v$index$is_outlier == 1L]
  expect_equal(v$curves[[i]]$phi_total, 1, tolerance = 1e-6)
  expect_false(isTRUE(all.equal(v$curves[[o]]$phi_total, 1, tolerance = 1e-3)))
})

test_that("windows are contiguous, in range, and capture essentially all mass", {
  for (k in seq_len(min(50, length(v$curves)))) {
    cv <- v$curves[[k]]
    expect_gte(cv$t0, CAL_GRID$t0)
    expect_lte(cv$t0 + (cv$len - 1) * CAL_GRID$dt, max(CAL_GRID$t))
    expect_equal(length(cv$L), cv$len)
    expect_equal(length(cv$Phi), cv$len)
    expect_equal(cv$Phi[1], 0)
    expect_equal(cv$Phi[cv$len], cv$phi_total, tolerance = 1e-12)
    expect_true(all(diff(cv$Phi) >= -1e-15))
  }
})

test_that("stan data has consistent ragged indexing", {
  pc <- site_phase_counts(d$dates)
  sd <- build_stan_data(d$dates, v, pc)
  expect_equal(sd$N, 164L)
  expect_equal(sd$S, 34L)
  expect_equal(sd$P, 74L)
  expect_equal(sum(sd$n_phase), sd$P)
  expect_equal(sd$phase_start[1], 1L)
  expect_equal(sd$phase_start + sd$n_phase - 1L, cumsum(sd$n_phase))
  expect_equal(sd$var_start + sd$n_var - 1L, cumsum(sd$n_var))
  expect_equal(length(sd$L_flat), sd$Ltot)
  expect_equal(sd$curve_pos + sd$curve_len - 1L, cumsum(sd$curve_len))
  expect_equal(sd$n_site_J2, 13L)
  expect_equal(sd$n_site_J4, 9L)
  expect_true(all(sd$class_id %in% 1:3))
})

test_that("indep and single configurations collapse phases correctly", {
  pc <- site_phase_counts(d$dates)
  si <- build_stan_data(d$dates, v, pc, model = "indep")
  expect_equal(si$P, 34L); expect_true(all(si$n_phase == 1L))
  ss <- build_stan_data(d$dates, v, pc, model = "single")
  expect_equal(ss$S, 1L); expect_equal(ss$P, 1L)
  expect_true(all(ss$site == 1L))
})

test_that("empty simplex-slot arrays stay typed integer for JSON serialisation", {
  # indep, single and sequential all leave J=1 everywhere, so both slot index
  # arrays are empty. A zero-length numeric or logical serialises to JSON in a
  # way CmdStan rejects for an int array, and the failure would surface only
  # after several model fits, so pin the type here.
  pc <- site_phase_counts(d$dates)
  for (m in c("indep", "single")) {
    sd <- build_stan_data(d$dates, v, pc, model = m)
    expect_equal(sd$n_site_J2, 0L)
    expect_equal(sd$n_site_J4, 0L)
    expect_length(sd$sites_J2, 0L)
    expect_length(sd$sites_J4, 0L)
    expect_type(sd$sites_J2, "integer")
    expect_type(sd$sites_J4, "integer")
    expect_type(sd$n_site_J2, "integer")
  }
})
