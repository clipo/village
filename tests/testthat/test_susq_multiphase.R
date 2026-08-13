source("../../R/susq_data.R")
source("../../R/susq_calibration.R")
source("../../R/susq_stan_data.R")
source("../../R/susq_fit.R")
suppressPackageStartupMessages(library(posterior))

MODEL_PHASES <- "../../stan/susq_phases.stan"
MODEL_SEQ <- "../../stan/susq_sequential.stan"

#' Simulate determinations from one or more known occupation windows at a site.
#' Each window is c(start_calBP, end_calBP) with start older than end.
simulate_site <- function(name, windows, n_each, err = 30, seed = 1) {
  set.seed(seed)
  cc <- intcal20_on_grid()
  rows <- list()
  for (w in seq_along(windows)) {
    tt <- runif(n_each[w], windows[[w]][2], windows[[w]][1])
    mu <- approx(cc$calBP, cc$c14, xout = tt)$y
    ce <- approx(cc$calBP, cc$err, xout = tt)$y
    rows[[w]] <- data.frame(
      site = name,
      c14_age = round(rnorm(length(tt), mu, sqrt(err^2 + ce^2))),
      c14_error = err, stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  n <- nrow(out)
  data.frame(
    det_id = seq_len(n), site = out$site, context = as.character(seq_len(n)),
    lab_no = paste0(name, "-", seq_len(n)),
    c14_age = out$c14_age, c14_error = out$c14_error,
    material_class = "short_lived", class_id = 1L, reference = "sim",
    pair_lab_no = NA_character_, pair_c14_age = NA_real_,
    pair_c14_error = NA_real_, pair_offset = NA_real_, preprocess = "none",
    stringsAsFactors = FALSE)
}

test_that("two well-separated episodes at one site are both recovered", {
  sim <- simulate_site("SIM", list(c(1000, 950), c(500, 450)), c(8, 8),
                       seed = 20260813)
  cc <- intcal20_on_grid()
  v <- build_variants(sim, cc)
  sd <- build_stan_data(sim, v, data.frame(site = "SIM", n = 16L, J = 4L))
  fit <- fit_susq(sd, dates = sim, model_file = MODEL_PHASES, refresh = 0)

  ident <- summarise_draws(
    fit$draws(c("phase_start_calBP", "phase_end_calBP")), "rhat")
  expect_lt(max(ident$rhat), 1.05)

  act <- fit$draws("active", format = "matrix")
  n_active <- rowSums(act)
  # The truth is two episodes. The posterior mode should be 2, and the model
  # must not routinely collapse to one.
  expect_equal(as.integer(names(which.max(table(n_active)))), 2L)
  expect_lt(mean(n_active == 1), 0.15)

  # Both true windows should be covered by an active phase in most draws.
  st <- fit$draws("phase_start_calBP", format = "matrix")
  en <- fit$draws("phase_end_calBP", format = "matrix")
  covers <- function(t) mean(rowSums((st >= t) & (en <= t) & (act == 1)) > 0)
  expect_gt(covers(975), 0.8)
  expect_gt(covers(475), 0.8)
  # And the long gap between them should mostly NOT be claimed as occupied.
  expect_lt(covers(725), 0.5)
})

test_that("a single episode is not split into spurious extra episodes", {
  sim <- simulate_site("SIM", list(c(700, 650)), 12, seed = 7)
  cc <- intcal20_on_grid()
  v <- build_variants(sim, cc)
  sd <- build_stan_data(sim, v, data.frame(site = "SIM", n = 12L, J = 4L))
  fit <- fit_susq(sd, dates = sim, model_file = MODEL_PHASES, refresh = 0)
  n_active <- rowSums(fit$draws("active", format = "matrix"))
  # The overfitted-mixture prior should empty the components data cannot
  # support rather than splitting one occupation across all four.
  expect_gt(mean(n_active <= 2), 0.8)
})

test_that("the sequential model runs and enforces non-overlap", {
  a <- simulate_site("A", list(c(900, 860)), 6, seed = 1)
  b <- simulate_site("B", list(c(500, 460)), 6, seed = 2)
  b$site <- "B"
  sim <- rbind(a, b)
  sim$det_id <- seq_len(nrow(sim))
  sim$lab_no <- paste0(sim$site, "-", sim$det_id)
  sim$context <- as.character(sim$det_id)

  cc <- intcal20_on_grid()
  v <- build_variants(sim, cc)
  sd <- build_stan_data(sim, v, site_phase_counts(sim), model = "sequential")
  fit <- fit_susq(sd, dates = sim, model_file = MODEL_SEQ, refresh = 0)

  st <- fit$draws("phase_start_calBP", format = "matrix")
  en <- fit$draws("phase_end_calBP", format = "matrix")
  # Phase 1 is the youngest. Its start must never reach phase 2's end, in
  # every draw, which is what "strictly serial" means.
  expect_true(all(st[, 1] <= en[, 2] + 1e-8))
  expect_true(all(st > en))
  # Sites arrive youngest first, so B (younger) must be phase 1.
  expect_equal(sd$site_names[1], "B")
})
