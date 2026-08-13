# Susquehanna Valley Contemporaneity Analysis Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Estimate how many settlements in the upper Susquehanna valley were occupied simultaneously, and for how long, from 164 radiocarbon determinations across 34 sites, using a Bayesian multi-phase occupation model in Stan.

**Architecture:** Latent calendar dates and discrete phase membership are both marginalized analytically. For every determination, the IntCal20 calibration likelihood is precomputed on a 1-year calendar grid along with its exact cumulative integral, so the probability of a determination given a uniform occupation phase `[b, a]` reduces to `(Phi(a) - Phi(b)) / (a - b)`, evaluated by a continuously differentiable interpolant inside Stan. Inbuilt age in wood and outlying measurements enter as additional precomputed likelihood variants weighted by estimated mixture probabilities. Each site carries an overfitted Dirichlet mixture over up to four candidate occupation phases, so the number of episodes is inferred rather than assumed. Occupation intervals and phase-activity indicators come out of `generated quantities`; all derived curves are computed in R from the draws.

**Tech Stack:** R 4.4, CmdStan 2.38.0 via `cmdstanr`, `IntCal` for IntCal20, `rcarbon` as an independent calibration check, `posterior`, `loo`, `ggplot2`, `systemfonts`, `testthat`.

**Spec:** `docs/superpowers/specs/2026-08-13-susquehanna-contemporaneity-design.md`

## Global Constraints

- Calendar time is **cal BP throughout**. Larger values are older. A phase runs from `a` (start, older, larger) to `b` (end, younger, smaller) with `a > b`. Convert to AD only at the figure layer, as `AD = 1950 - calBP`.
- Calendar grid: 1-year spacing, `cal BP 0` to `cal BP 2000` inclusive (`G = 2001`).
- Inbuilt-age scales `kappa`: exactly `c(0, 15, 30, 60, 120, 240)` years, in that order. Index 1 is always `kappa = 0`.
- Outlier error inflation factor: exactly `5`.
- Material classes are exactly `short_lived`, `wood`, `indeterminate`, encoded as integers `1`, `2`, `3` in that order.
- Phase counts: `J_s = 4` if the site has >= 6 determinations, `2` if 3 to 5, `1` if 1 or 2. Totals 74 phases across 34 sites.
- Dirichlet concentration for phase weights: `alpha = 1`, so each site's simplex prior is `Dirichlet(1 / J_s)`.
- Duration prior: `dur ~ Lognormal(mu_d, sigma_d)`, `mu_d ~ Normal(log 60, 1)`, `sigma_d ~ Normal(0, 1) T[0,]`. Sensitivity runs replace `log 60` with `log 25` and `log 150`, and nothing else changes.
- Outlier rate prior: `rho_c ~ Beta(2, 18)` for each of the three classes.
- All likelihood variants of one determination share a single normalising constant. **Never normalise variants separately**; their relative scale carries the outlier and inbuilt-age evidence.
- Random seed `20260813` everywhere. Four chains, 2000 warmup, 2000 sampling, `adapt_delta = 0.95`, `max_treedepth = 12`.
- Figures: PNG via `ragg`, 300 dpi, 7 in wide (3.5 for single column), Arial via `systemfonts`, Okabe-Ito categorical, viridis continuous, no titles, minimal gridlines.
- Every script sources its dependencies explicitly and is runnable standalone from the project root.
- `renv::snapshot()` after any new package is added.
- Outputs go to `output/susquehanna/`. Figures go to `output/susquehanna/figures/`.
- Sourcing `PATH` for GDAL is required before anything loads `rcarbon`: `/opt/homebrew/opt/{gdal,geos,proj}/bin`.

---

### Task 1: Data preparation

**Files:**
- Create: `R/susq_data.R`
- Test: `tests/testthat/test_susq_data.R`
- Read only: `data/Local Dates.xlsx`, `data/site_parse.csv`, `data/material_classes.csv`

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `prepare_susq_data(xlsx = "data/Local Dates.xlsx", site_map = "data/site_parse.csv", material_map = "data/material_classes.csv")` returns a list with `$dates` (data frame, one row per modeled determination) and `$quality` (data frame of flags).
  - `$dates` columns, in this order: `det_id` (integer 1..N), `site` (character), `context` (character), `lab_no` (character), `c14_age` (numeric), `c14_error` (numeric), `material_class` (character, one of `short_lived`/`wood`/`indeterminate`), `class_id` (integer 1/2/3), `reference` (character), `pair_lab_no` (character or NA), `pair_c14_age` (numeric or NA), `pair_c14_error` (numeric or NA), `pair_offset` (numeric or NA, years the paired member is younger), `preprocess` (character: `none`, `pooled_replicate`, `wiggle_match`).
  - `site_phase_counts(dates)` returns a data frame with `site`, `n`, `J`, ordered by site name.

- [ ] **Step 1: Write the failing test**

```r
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test_susq_data.R")'`
Expected: FAIL, `could not find function "prepare_susq_data"`.

- [ ] **Step 3: Write the implementation**

```r
# R/susq_data.R
#' Load and preprocess the Susquehanna radiocarbon determinations.
#'
#' Applies the audited site and material lookup tables, pools the failed
#' Bates_66 replicate pair, and folds the two ring-counted wiggle-match pairs
#' into composite determinations. See spec sections 2.2 to 2.5.

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(stringr)
})

MATERIAL_CLASSES <- c("short_lived", "wood", "indeterminate")

#' Ward and Wilson (1978) test for a set of determinations on one sample.
#' Returns the pooled age, the pooled error, the test statistic and its df.
ward_wilson_pool <- function(ages, errors) {
  w <- 1 / errors^2
  pooled <- sum(ages * w) / sum(w)
  se <- sqrt(1 / sum(w))
  tstat <- sum(w * (ages - pooled)^2)
  df <- length(ages) - 1L
  crit <- stats::qchisq(0.95, df)
  # When the determinations are not statistically consistent, inflate the
  # pooled error by sqrt(T / df) rather than reporting a spuriously precise
  # combined age.
  if (tstat > crit) se <- se * sqrt(tstat / df)
  list(age = pooled, error = se, t = tstat, df = df,
       consistent = tstat <= crit)
}

prepare_susq_data <- function(xlsx = "data/Local Dates.xlsx",
                              site_map = "data/site_parse.csv",
                              material_map = "data/material_classes.csv") {
  raw <- readxl::read_excel(xlsx)
  names(raw)[1:6] <- c("site_raw", "material_raw", "lab_no",
                       "c14_age", "c14_error", "reference")
  raw <- as.data.frame(raw[, 1:6])
  raw$c14_age <- as.numeric(raw$c14_age)
  raw$c14_error <- as.numeric(raw$c14_error)

  smap <- read.csv(site_map, stringsAsFactors = FALSE)
  mmap <- read.csv(material_map, stringsAsFactors = FALSE)

  d <- raw %>%
    left_join(smap[, c("site_raw", "site", "context")], by = "site_raw") %>%
    left_join(mmap[, c("material_raw", "material_class", "ring_lo", "ring_hi")],
              by = "material_raw")

  stopifnot(!any(is.na(d$site)), !any(is.na(d$material_class)))

  quality <- collect_quality_flags(d)

  d <- pool_replicates(d)
  d <- fold_wiggle_matches(d)

  d$class_id <- match(d$material_class, MATERIAL_CLASSES)
  d <- d[order(d$site, d$context, d$lab_no), ]
  d$det_id <- seq_len(nrow(d))

  keep <- c("det_id", "site", "context", "lab_no", "c14_age", "c14_error",
            "material_class", "class_id", "reference",
            "pair_lab_no", "pair_c14_age", "pair_c14_error", "pair_offset",
            "preprocess")
  list(dates = d[, keep], quality = quality)
}

#' Pool replicate determinations. A lab code ending in "r" denotes a rerun of
#' the sample identified by the code without the suffix.
pool_replicates <- function(d) {
  d$preprocess <- "none"
  reps <- grep("[0-9]r$", d$lab_no, value = TRUE)
  for (r in reps) {
    base <- sub("r$", "", r)
    idx <- which(d$lab_no %in% c(base, r))
    if (length(idx) < 2) next
    p <- ward_wilson_pool(d$c14_age[idx], d$c14_error[idx])
    keep <- idx[d$lab_no[idx] == base]
    drop <- setdiff(idx, keep)
    d$c14_age[keep] <- p$age
    d$c14_error[keep] <- p$error
    d$preprocess[keep] <- "pooled_replicate"
    d <- d[-drop, , drop = FALSE]
  }
  d
}

#' Fold ring-counted pairs into one composite determination. The member with
#' the lower ring numbers is retained as the reference; the other is recorded
#' as the pair, younger by the difference of the ring-span midpoints.
fold_wiggle_matches <- function(d) {
  d$pair_lab_no <- NA_character_
  d$pair_c14_age <- NA_real_
  d$pair_c14_error <- NA_real_
  d$pair_offset <- NA_real_

  keyed <- which(!is.na(d$ring_lo))
  if (!length(keyed)) return(d)
  groups <- split(keyed, paste(d$site[keyed], d$context[keyed]))
  drop <- integer(0)
  for (g in groups) {
    if (length(g) != 2L) next
    mid <- (d$ring_lo[g] + d$ring_hi[g]) / 2
    ref <- g[which.min(mid)]
    oth <- g[which.max(mid)]
    d$pair_lab_no[ref] <- d$lab_no[oth]
    d$pair_c14_age[ref] <- d$c14_age[oth]
    d$pair_c14_error[ref] <- d$c14_error[oth]
    d$pair_offset[ref] <- max(mid) - min(mid)
    d$preprocess[ref] <- "wiggle_match"
    drop <- c(drop, oth)
  }
  if (length(drop)) d <- d[-drop, , drop = FALSE]
  d
}

collect_quality_flags <- function(d) {
  flags <- list()
  add <- function(kind, detail) flags[[length(flags) + 1L]] <<-
    data.frame(kind = kind, detail = detail, stringsAsFactors = FALSE)

  dup <- d$lab_no[duplicated(d$lab_no) & !is.na(d$lab_no)]
  for (l in unique(dup)) {
    rows <- d[d$lab_no == l, ]
    add("duplicate_lab_code",
        sprintf("%s appears %d times: %s",
                l, nrow(rows),
                paste(sprintf("%s %g+/-%g", rows$site, rows$c14_age, rows$c14_error),
                      collapse = "; ")))
  }
  for (i in which(grepl("none provided|\\?", d$lab_no))) {
    add("questionable_lab_code",
        sprintf("%s at %s (%g+/-%g)", d$lab_no[i], d$site[i],
                d$c14_age[i], d$c14_error[i]))
  }
  do.call(rbind, flags)
}

site_phase_counts <- function(dates) {
  pc <- as.data.frame(table(dates$site), stringsAsFactors = FALSE)
  names(pc) <- c("site", "n")
  pc$J <- ifelse(pc$n >= 6L, 4L, ifelse(pc$n >= 3L, 2L, 1L))
  pc[order(pc$site), ]
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test_susq_data.R")'`
Expected: PASS, 6 tests.

If the material-class counts differ from the test, do **not** change the test to match. Re-read `data/material_classes.csv`, confirm the audited assignment, and fix the code.

- [ ] **Step 5: Commit**

```bash
git add R/susq_data.R tests/testthat/test_susq_data.R
git commit -m "Add data preparation with replicate pooling and wiggle-match folding"
```

---

### Task 2: Calibration likelihood and its cumulative integral

**Files:**
- Create: `R/susq_calibration.R`
- Test: `tests/testthat/test_susq_calibration.R`

**Interfaces:**
- Consumes: `prepare_susq_data()` from Task 1.
- Produces:
  - `CAL_GRID` constant: `list(t0 = 0, dt = 1, G = 2001, t = 0:2000)`.
  - `intcal20_on_grid()` returns a data frame with `calBP`, `c14`, `err`, interpolated to `CAL_GRID$t`.
  - `calib_likelihood(c14_age, c14_error, curve, inflate = 1)` returns a numeric vector of length `G`, the unnormalised `p(y | t)`.
  - `convolve_inbuilt(L, kappa, dt)` returns a vector of length `G`, `L` convolved with an exponential inbuilt-age distribution of mean `kappa`. `kappa = 0` returns `L` unchanged.
  - `pair_likelihood(L_ref, L_other, offset, dt)` returns the elementwise product of `L_ref` with `L_other` shifted so that the paired member sits `offset` years younger.
  - `cumulative_integral(L, dt)` returns a vector of length `G` giving the exact integral of the piecewise-linear interpolant of `L` from the grid start to each knot. Element 1 is 0.

- [ ] **Step 1: Write the failing test**

```r
# tests/testthat/test_susq_calibration.R
source("../../R/susq_calibration.R")

test_that("the calendar grid matches the global constraint", {
  expect_equal(CAL_GRID$t0, 0); expect_equal(CAL_GRID$dt, 1)
  expect_equal(CAL_GRID$G, 2001L); expect_equal(length(CAL_GRID$t), 2001L)
})

test_that("IntCal20 is interpolated onto the grid without gaps", {
  cc <- intcal20_on_grid()
  expect_equal(nrow(cc), 2001L)
  expect_false(any(is.na(cc$c14))); expect_false(any(is.na(cc$err)))
  expect_true(all(cc$err > 0))
  expect_equal(cc$calBP, CAL_GRID$t)
})

test_that("cumulative_integral integrates a piecewise-linear function exactly", {
  # L(t) = t on a unit grid: integral from 0 to k is k^2 / 2, exactly
  # reproduced by the trapezoid rule on a linear function.
  L <- as.numeric(0:10)
  got <- cumulative_integral(L, dt = 1)
  expect_equal(got, (0:10)^2 / 2, tolerance = 1e-12)
  expect_equal(got[1], 0)
})

test_that("cumulative_integral is monotone non-decreasing for a likelihood", {
  cc <- intcal20_on_grid()
  L <- calib_likelihood(700, 30, cc)
  Phi <- cumulative_integral(L, CAL_GRID$dt)
  expect_true(all(diff(Phi) >= -1e-15))
})

test_that("calibrated density agrees with rcarbon", {
  # rcarbon is an independent implementation. Compare the normalised density
  # over the grid; agreement to 1e-3 in total variation is ample.
  cc <- intcal20_on_grid()
  L <- calib_likelihood(700, 30, cc)
  ours <- L / sum(L)

  rc <- rcarbon::calibrate(700, 30, calCurves = "intcal20",
                           normalised = TRUE, verbose = FALSE)
  g <- rc$grids[[1]]
  theirs <- rep(0, CAL_GRID$G)
  keep <- g$calBP >= 0 & g$calBP <= 2000
  theirs[match(g$calBP[keep], CAL_GRID$t)] <- g$PrDens[keep]
  theirs <- theirs / sum(theirs)

  expect_lt(0.5 * sum(abs(ours - theirs)), 1e-3)
})

test_that("exponential convolution shifts mass older and preserves total mass", {
  cc <- intcal20_on_grid()
  L <- calib_likelihood(700, 30, cc)
  Lk <- convolve_inbuilt(L, kappa = 60, dt = CAL_GRID$dt)
  expect_equal(convolve_inbuilt(L, 0, CAL_GRID$dt), L)
  expect_equal(sum(Lk), sum(L), tolerance = 1e-6 * sum(L))
  # Inbuilt age means the sample grew before deposition, so the event-date
  # likelihood moves toward younger calendar dates, i.e. smaller cal BP.
  centroid <- function(v) sum(v * CAL_GRID$t) / sum(v)
  expect_lt(centroid(Lk), centroid(L))
  expect_gt(centroid(L) - centroid(Lk), 30)
  expect_lt(centroid(L) - centroid(Lk), 90)
})

test_that("pair_likelihood sharpens a wiggle-matched pair", {
  cc <- intcal20_on_grid()
  L1 <- calib_likelihood(870, 25, cc)
  L2 <- calib_likelihood(820, 25, cc)
  Lp <- pair_likelihood(L1, L2, offset = 17, dt = CAL_GRID$dt)
  sd_of <- function(v) { p <- v / sum(v); m <- sum(p * CAL_GRID$t)
                         sqrt(sum(p * (CAL_GRID$t - m)^2)) }
  expect_lt(sd_of(Lp), sd_of(L1))
  expect_true(all(Lp >= 0))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test_susq_calibration.R")'`
Expected: FAIL, `object 'CAL_GRID' not found`.

- [ ] **Step 3: Write the implementation**

```r
# R/susq_calibration.R
#' Calibration machinery. Builds, for each determination, the likelihood of the
#' observed radiocarbon age as a function of calendar date, on a fixed 1-year
#' grid, together with the exact integral of its piecewise-linear interpolant.
#' See spec section 3.3.

CAL_GRID <- list(t0 = 0, dt = 1, G = 2001L, t = 0:2000)

#' IntCal20 interpolated onto the analysis grid.
intcal20_on_grid <- function() {
  cc <- IntCal::ccurve(cc = 1)
  colnames(cc) <- c("calBP", "c14", "err")
  cc <- as.data.frame(cc)
  cc <- cc[cc$calBP <= max(CAL_GRID$t) + 50, ]
  data.frame(
    calBP = CAL_GRID$t,
    c14   = stats::approx(cc$calBP, cc$c14, xout = CAL_GRID$t, rule = 2)$y,
    err   = stats::approx(cc$calBP, cc$err, xout = CAL_GRID$t, rule = 2)$y
  )
}

#' Unnormalised p(y | t): the probability of the observed radiocarbon age given
#' a true calendar date, combining measurement and curve error.
#' `inflate` multiplies the measurement error, for the outlier variant.
calib_likelihood <- function(c14_age, c14_error, curve, inflate = 1) {
  s <- sqrt((c14_error * inflate)^2 + curve$err^2)
  stats::dnorm(c14_age, mean = curve$c14, sd = s)
}

#' Convolve a likelihood with an exponential inbuilt-age distribution of mean
#' `kappa`. A sample with inbuilt age d grew d years before deposition, so the
#' event-date likelihood at t is the likelihood of the growth date at t + d
#' averaged over d. Mass therefore moves toward smaller cal BP.
convolve_inbuilt <- function(L, kappa, dt) {
  if (kappa <= 0) return(L)
  # Kernel over d = 0, dt, 2dt, ... truncated where it contributes nothing.
  dmax <- ceiling(qexp(1 - 1e-8, rate = 1 / kappa) / dt) * dt
  d <- seq(0, dmax, by = dt)
  w <- stats::dexp(d, rate = 1 / kappa)
  w <- w / sum(w)
  G <- length(L)
  out <- numeric(G)
  for (k in seq_along(d)) {
    shift <- (k - 1L)
    if (shift == 0L) {
      out <- out + w[k] * L
    } else {
      # event date t maps to growth date t + shift, so read L from higher index
      idx <- seq_len(G - shift)
      out[idx] <- out[idx] + w[k] * L[idx + shift]
    }
  }
  out
}

#' Likelihood for a wiggle-matched pair. `offset` is how many years younger the
#' paired member is than the reference member, known from ring counting. The
#' returned curve is indexed by the calendar date of the reference member.
pair_likelihood <- function(L_ref, L_other, offset, dt) {
  shift <- as.integer(round(offset / dt))
  G <- length(L_ref)
  shifted <- numeric(G)
  # The paired member's date is (reference date - offset), i.e. a lower cal BP,
  # so its likelihood is read at index i - shift.
  idx <- seq.int(shift + 1L, G)
  shifted[idx] <- L_other[idx - shift]
  L_ref * shifted
}

#' Exact integral of the piecewise-linear interpolant of L, from the grid start
#' to each knot. Element 1 is 0. This is the trapezoid rule, which is exact for
#' a piecewise-linear integrand.
cumulative_integral <- function(L, dt) {
  c(0, cumsum((utils::head(L, -1) + utils::tail(L, -1)) / 2 * dt))
}
```

- [ ] **Step 4: Run test to verify it passes**

Run:
```bash
export PATH="/opt/homebrew/opt/gdal/bin:/opt/homebrew/opt/geos/bin:/opt/homebrew/opt/proj/bin:$PATH"
Rscript -e 'testthat::test_file("tests/testthat/test_susq_calibration.R")'
```
Expected: PASS, 7 tests. The rcarbon comparison is the one that matters most; if it fails, the calibration is wrong and nothing downstream is trustworthy.

- [ ] **Step 5: Commit**

```bash
git add R/susq_calibration.R tests/testthat/test_susq_calibration.R
git commit -m "Add IntCal20 calibration likelihood, inbuilt-age convolution and cumulative integral"
```

---

### Task 3: Pack likelihood variants into ragged Stan data

**Files:**
- Create: `R/susq_stan_data.R`
- Test: `tests/testthat/test_susq_stan_data.R`

**Interfaces:**
- Consumes: `prepare_susq_data()` (Task 1), everything in `R/susq_calibration.R` (Task 2).
- Produces:
  - `build_variants(dates, curve)` returns a list with `curves` (list of per-variant lists each holding `L`, `Phi`, `t0`, `len`, `phi_total`) and `index` (data frame with `variant_id`, `det_id`, `kappa_idx` (0 for the outlier variant), `is_outlier`).
  - `build_stan_data(dates, variants, phase_counts, mu_d_prior = log(60), model = c("multi", "indep", "single"))` returns the named list passed to CmdStan.
  - `KAPPAS` constant `c(0, 15, 30, 60, 120, 240)`.
  - `OUTLIER_INFLATE` constant `5`.

Variants are windowed: only the stretch of grid where the curve carries meaningful mass is stored, which keeps the Stan data around 5 MB rather than 25 MB. Outside its window a variant's cumulative integral is 0 below and `phi_total` above.

- [ ] **Step 1: Write the failing test**

```r
# tests/testthat/test_susq_stan_data.R
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
  expect_equal(sd$phase_start + sd$n_phase - 1L,
               cumsum(sd$n_phase))
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test_susq_stan_data.R")'`
Expected: FAIL, `could not find function "build_variants"`.

- [ ] **Step 3: Write the implementation**

```r
# R/susq_stan_data.R
#' Build the ragged likelihood-variant structure and assemble Stan data.
#' See spec sections 3.3 to 3.6.

KAPPAS <- c(0, 15, 30, 60, 120, 240)
OUTLIER_INFLATE <- 5
WINDOW_REL_TOL <- 1e-10
WINDOW_PAD <- 25L

#' Trim a curve to the stretch carrying meaningful mass, and return the window
#' along with the exact cumulative integral inside it.
window_curve <- function(L, dt, t_grid) {
  thresh <- max(L) * WINDOW_REL_TOL
  keep <- which(L > thresh)
  lo <- max(1L, min(keep) - WINDOW_PAD)
  hi <- min(length(L), max(keep) + WINDOW_PAD)
  Lw <- L[lo:hi]
  Phi <- cumulative_integral(Lw, dt)
  list(L = Lw, Phi = Phi, t0 = t_grid[lo], len = length(Lw),
       phi_total = Phi[length(Phi)])
}

build_variants <- function(dates, curve) {
  curves <- list(); index <- list(); vid <- 0L

  for (i in seq_len(nrow(dates))) {
    row <- dates[i, ]
    base <- calib_likelihood(row$c14_age, row$c14_error, curve)
    out  <- calib_likelihood(row$c14_age, row$c14_error, curve,
                             inflate = OUTLIER_INFLATE)

    # A wiggle-matched pair contributes the product of its members' likelihoods
    # with the known ring offset applied, in both the base and outlier forms.
    if (!is.na(row$pair_offset)) {
      p_base <- calib_likelihood(row$pair_c14_age, row$pair_c14_error, curve)
      p_out  <- calib_likelihood(row$pair_c14_age, row$pair_c14_error, curve,
                                 inflate = OUTLIER_INFLATE)
      base <- pair_likelihood(base, p_base, row$pair_offset, CAL_GRID$dt)
      out  <- pair_likelihood(out,  p_out,  row$pair_offset, CAL_GRID$dt)
    }

    # One normaliser for every variant of this determination, so that the
    # relative scale between variants is preserved. See spec section 3.4.
    norm <- sum((utils::head(base, -1) + utils::tail(base, -1)) / 2 * CAL_GRID$dt)
    stopifnot(norm > 0)
    base <- base / norm
    out  <- out  / norm

    n_kappa <- if (row$material_class == "short_lived") 1L else length(KAPPAS)
    for (k in seq_len(n_kappa)) {
      Lk <- convolve_inbuilt(base, KAPPAS[k], CAL_GRID$dt)
      vid <- vid + 1L
      curves[[vid]] <- window_curve(Lk, CAL_GRID$dt, CAL_GRID$t)
      index[[vid]] <- data.frame(variant_id = vid, det_id = row$det_id,
                                 kappa_idx = k, is_outlier = 0L)
    }
    vid <- vid + 1L
    curves[[vid]] <- window_curve(out, CAL_GRID$dt, CAL_GRID$t)
    index[[vid]] <- data.frame(variant_id = vid, det_id = row$det_id,
                               kappa_idx = 0L, is_outlier = 1L)
  }
  list(curves = curves, index = do.call(rbind, index))
}

build_stan_data <- function(dates, variants, phase_counts,
                            mu_d_prior = log(60),
                            model = c("multi", "indep", "single")) {
  model <- match.arg(model)

  sites <- sort(unique(dates$site))
  site_id <- match(dates$site, sites)
  J <- phase_counts$J[match(sites, phase_counts$site)]

  if (model == "indep") J <- rep(1L, length(sites))
  if (model == "single") { sites <- "ALL"; site_id <- rep(1L, nrow(dates)); J <- 1L }

  S <- length(J); P <- sum(J)
  phase_start <- as.integer(cumsum(c(1L, utils::head(J, -1))))

  idx <- variants$index[order(variants$index$det_id, variants$index$variant_id), ]
  n_var <- as.integer(table(factor(idx$det_id, levels = dates$det_id)))
  var_start <- as.integer(cumsum(c(1L, utils::head(n_var, -1))))

  ord <- idx$variant_id
  cvs <- variants$curves[ord]
  curve_len <- vapply(cvs, function(x) as.integer(x$len), 1L)
  curve_pos <- as.integer(cumsum(c(1L, utils::head(curve_len, -1))))

  list(
    N = nrow(dates), S = S, P = P,
    site = site_id,
    phase_start = phase_start, n_phase = as.integer(J),
    class_id = as.integer(dates$class_id),

    V = nrow(idx),
    var_start = var_start, n_var = n_var,
    kappa_idx = as.integer(idx$kappa_idx),
    is_outlier = as.integer(idx$is_outlier),

    Ltot = sum(curve_len),
    L_flat = unlist(lapply(cvs, `[[`, "L"), use.names = FALSE),
    Phi_flat = unlist(lapply(cvs, `[[`, "Phi"), use.names = FALSE),
    curve_pos = curve_pos, curve_len = curve_len,
    curve_t0 = vapply(cvs, function(x) as.numeric(x$t0), 0),
    phi_total = vapply(cvs, function(x) as.numeric(x$phi_total), 0),
    dt = CAL_GRID$dt,

    n_kappa = length(KAPPAS),
    cal_min = min(CAL_GRID$t), cal_max = max(CAL_GRID$t),

    sites_J2 = which(J == 2L), n_site_J2 = sum(J == 2L),
    sites_J4 = which(J == 4L), n_site_J4 = sum(J == 4L),

    mu_d_prior_mean = mu_d_prior,
    site_names = sites
  )
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test_susq_stan_data.R")'`
Expected: PASS, 6 tests.

- [ ] **Step 5: Commit**

```bash
git add R/susq_stan_data.R tests/testthat/test_susq_stan_data.R
git commit -m "Add ragged likelihood-variant packing and Stan data assembly"
```

---

### Task 4: The interpolant, tested inside Stan

**Files:**
- Create: `stan/susq_functions.stan`
- Create: `stan/test_phi_at.stan`
- Test: `tests/testthat/test_phi_at.R`

**Interfaces:**
- Produces: Stan function `phi_at(real t, int v, vector L_flat, vector Phi_flat, array[] int curve_pos, array[] int curve_len, array[] real curve_t0, array[] real phi_total, real dt)`, returning the cumulative integral of variant `v` evaluated at calendar date `t`, clamped to `0` below the window and `phi_total` above it. Included by both Stan models via `#include susq_functions.stan`.

This task exists on its own because the interpolant is the numerical heart of the model. If it is wrong, every posterior is wrong and no other test will say so clearly.

- [ ] **Step 1: Write the Stan function file and the test harness model**

```stan
// stan/susq_functions.stan
/**
 * Cumulative integral of a variant's likelihood curve, evaluated at an
 * arbitrary calendar date.
 *
 * Each variant stores its likelihood L on a uniform grid of `curve_len[v]`
 * knots starting at cal BP `curve_t0[v]`, together with the exact integral
 * Phi of the piecewise-linear interpolant of L at each knot. Because L is
 * treated as piecewise linear, Phi is piecewise quadratic and therefore
 * continuously differentiable in t, which is what Hamiltonian Monte Carlo
 * needs. Outside the stored window the curve carries no meaningful mass, so
 * Phi clamps to 0 below and to phi_total above.
 */
real phi_at(real t, int v,
            vector L_flat, vector Phi_flat,
            array[] int curve_pos, array[] int curve_len,
            array[] real curve_t0, array[] real phi_total,
            real dt) {
  real lo = curve_t0[v];
  real hi = curve_t0[v] + (curve_len[v] - 1) * dt;
  if (t <= lo) return 0;
  if (t >= hi) return phi_total[v];
  real x = (t - lo) / dt;
  int g = 1 + to_int(floor(x));
  real h = t - (lo + (g - 1) * dt);
  int p = curve_pos[v] + g - 1;
  return Phi_flat[p] + L_flat[p] * h
         + (L_flat[p + 1] - L_flat[p]) * h * h / (2 * dt);
}

/**
 * Probability of a determination given a uniform occupation phase [b, a],
 * with the latent calendar date integrated out. `a` is the older boundary.
 */
real phase_average(real a, real b, int v,
                   vector L_flat, vector Phi_flat,
                   array[] int curve_pos, array[] int curve_len,
                   array[] real curve_t0, array[] real phi_total,
                   real dt) {
  return (phi_at(a, v, L_flat, Phi_flat, curve_pos, curve_len,
                 curve_t0, phi_total, dt)
        - phi_at(b, v, L_flat, Phi_flat, curve_pos, curve_len,
                 curve_t0, phi_total, dt)) / (a - b);
}
```

```stan
// stan/test_phi_at.stan
// Harness that exposes the interpolant to R for direct numerical comparison.
functions {
#include susq_functions.stan
}
data {
  int<lower=1> Ltot;
  vector[Ltot] L_flat;
  vector[Ltot] Phi_flat;
  int<lower=1> V;
  array[V] int curve_pos;
  array[V] int curve_len;
  array[V] real curve_t0;
  array[V] real phi_total;
  real dt;
  int<lower=1> M;
  array[M] real t_query;
  array[M] int v_query;
}
parameters { real dummy; }
model { dummy ~ std_normal(); }
generated quantities {
  array[M] real phi_out;
  for (m in 1:M)
    phi_out[m] = phi_at(t_query[m], v_query[m], L_flat, Phi_flat,
                        curve_pos, curve_len, curve_t0, phi_total, dt);
}
```

- [ ] **Step 2: Write the failing test**

```r
# tests/testthat/test_phi_at.R
source("../../R/susq_data.R")
source("../../R/susq_calibration.R")
source("../../R/susq_stan_data.R")
library(cmdstanr)

test_that("the Stan interpolant matches R integration of the same curve", {
  d <- prepare_susq_data(); cc <- intcal20_on_grid()
  v <- build_variants(d$dates[1:5, ], cc)
  sd <- build_stan_data(d$dates[1:5, ], v,
                        site_phase_counts(d$dates[1:5, ]))

  # Query each variant at 200 dates spanning and overrunning its window,
  # so the clamping branches are exercised too.
  set.seed(20260813)
  M <- 200L
  vq <- sample(seq_len(sd$V), M, replace = TRUE)
  lo <- sd$curve_t0[vq]
  hi <- lo + (sd$curve_len[vq] - 1) * sd$dt
  tq <- runif(M, lo - 40, hi + 40)

  mod <- cmdstan_model("stan/test_phi_at.stan", include_paths = "stan")
  fit <- mod$sample(data = c(sd[c("Ltot","L_flat","Phi_flat","curve_pos",
                                  "curve_len","curve_t0","phi_total","dt")],
                             list(V = sd$V, M = M, t_query = tq, v_query = vq)),
                    chains = 1, iter_warmup = 100, iter_sampling = 1,
                    seed = 20260813, refresh = 0)
  got <- as.numeric(fit$draws("phi_out", format = "matrix")[1, ])

  # Reference: integrate the stored piecewise-linear curve in R.
  ref <- vapply(seq_len(M), function(m) {
    vv <- vq[m]; p0 <- sd$curve_pos[vv]; n <- sd$curve_len[vv]
    L <- sd$L_flat[p0:(p0 + n - 1)]; Phi <- sd$Phi_flat[p0:(p0 + n - 1)]
    t0 <- sd$curve_t0[vv]; t <- tq[m]
    if (t <= t0) return(0)
    if (t >= t0 + (n - 1) * sd$dt) return(sd$phi_total[vv])
    g <- 1L + as.integer(floor((t - t0) / sd$dt))
    h <- t - (t0 + (g - 1) * sd$dt)
    Phi[g] + L[g] * h + (L[g + 1] - L[g]) * h^2 / (2 * sd$dt)
  }, 0)

  expect_equal(got, ref, tolerance = 1e-9)
})

test_that("the interpolant reproduces fine-grained numerical integration", {
  # Independent check that does not reuse the same formula: integrate the
  # underlying likelihood on a 0.01-year grid and compare.
  cc <- intcal20_on_grid()
  L <- calib_likelihood(700, 30, cc)
  L <- L / sum((head(L, -1) + tail(L, -1)) / 2)
  Phi <- cumulative_integral(L, 1)
  for (t in c(120.3, 450.7, 661.25, 900.5)) {
    fine <- seq(0, t, by = 0.01)
    dense <- approx(CAL_GRID$t, L, xout = fine)$y
    num <- sum((head(dense, -1) + tail(dense, -1)) / 2 * 0.01)
    g <- 1L + as.integer(floor(t)); h <- t - (g - 1)
    ours <- Phi[g] + L[g] * h + (L[g + 1] - L[g]) * h^2 / 2
    expect_equal(ours, num, tolerance = 1e-6)
  }
})
```

- [ ] **Step 3: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test_phi_at.R")'`
Expected: FAIL, the model file does not compile or does not exist.

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test_phi_at.R")'`
Expected: PASS, 2 tests. Compilation takes about a minute the first time.

- [ ] **Step 5: Commit**

```bash
git add stan/susq_functions.stan stan/test_phi_at.stan tests/testthat/test_phi_at.R
git commit -m "Add Stan phase-average interpolant with direct numerical verification"
```

---

### Task 5: The multi-phase occupation model

**Files:**
- Create: `stan/susq_phases.stan`
- Create: `R/susq_fit.R`
- Test: `tests/testthat/test_susq_phases_single.R`

**Interfaces:**
- Consumes: `stan/susq_functions.stan` (Task 4), `build_stan_data()` (Task 3).
- Produces:
  - `fit_susq(stan_data, model_file = "stan/susq_phases.stan", ...)` returns a `CmdStanMCMC` object, with initialisation drawn from the site's own calibrated dates.
  - `init_from_data(stan_data, dates, chains)` returns a list of init lists.
  - `check_diagnostics(fit)` returns a data frame with `max_rhat`, `min_ess_bulk`, `min_ess_tail`, `n_divergent`, `ebfmi_min`, and stops with an informative error if the spec's thresholds are breached and `strict = TRUE`.

- [ ] **Step 1: Write the Stan model**

```stan
// stan/susq_phases.stan
/**
 * Multi-phase occupation model for radiocarbon determinations grouped by site.
 *
 * Each site carries up to four candidate occupation phases. Determinations
 * are uniformly distributed within whichever phase they belong to, and phase
 * membership is marginalised as a Dirichlet mixture. The latent calendar date
 * of each determination is marginalised analytically using precomputed
 * cumulative calibration integrals, so no per-determination date parameter
 * appears. Inbuilt age and measurement outliers enter as a second mixture over
 * precomputed likelihood variants.
 *
 * Setting every n_phase to 1 gives the independent-phase model; collapsing all
 * determinations onto one site gives the single-phase model.
 */
functions {
#include susq_functions.stan
}

data {
  int<lower=1> N;                                  // determinations
  int<lower=1> S;                                  // sites
  int<lower=1> P;                                  // candidate phases in total
  array[N] int<lower=1, upper=S> site;
  array[S] int<lower=1, upper=P> phase_start;
  array[S] int<lower=1> n_phase;
  array[N] int<lower=1, upper=3> class_id;         // 1 short_lived, 2 wood, 3 indeterminate

  int<lower=1> V;                                  // likelihood variants in total
  array[N] int<lower=1, upper=V> var_start;
  array[N] int<lower=1> n_var;
  array[V] int<lower=0> kappa_idx;                 // 0 marks the outlier variant
  array[V] int<lower=0, upper=1> is_outlier;

  int<lower=1> Ltot;
  vector[Ltot] L_flat;
  vector[Ltot] Phi_flat;
  array[V] int<lower=1> curve_pos;
  array[V] int<lower=1> curve_len;
  array[V] real curve_t0;
  array[V] real phi_total;
  real<lower=0> dt;

  int<lower=1> n_kappa;
  real cal_min;
  real cal_max;

  int<lower=0> n_site_J2;
  array[n_site_J2] int<lower=1, upper=S> sites_J2;
  int<lower=0> n_site_J4;
  array[n_site_J4] int<lower=1, upper=S> sites_J4;

  real mu_d_prior_mean;
  int<lower=0, upper=1> prior_only;                // 1 runs the prior predictive
}

transformed data {
  // Map each site to its slot in the J=2 or J=4 simplex arrays. Sites with a
  // single phase need no simplex.
  array[S] int slot2 = rep_array(0, S);
  array[S] int slot4 = rep_array(0, S);
  for (m in 1:n_site_J2) slot2[sites_J2[m]] = m;
  for (m in 1:n_site_J4) slot4[sites_J4[m]] = m;
}

parameters {
  // Oldest phase midpoint at each site, then positive gaps to the younger ones.
  vector<lower=cal_min, upper=cal_max>[S] mid_first;
  vector[P - S] log_gap;

  // Durations, non-centred on a hierarchical lognormal.
  vector[P] log_dur_raw;
  real mu_d;
  real<lower=0> sigma_d;

  // Phase weights.
  array[n_site_J2] simplex[2] pi2;
  array[n_site_J4] simplex[4] pi4;

  // Inbuilt-age distributions for wood and indeterminate material.
  array[2] simplex[n_kappa] q;

  // Outlier rate per material class.
  vector<lower=0, upper=1>[3] rho;
}

transformed parameters {
  vector[P] mid;
  vector[P] dur = exp(mu_d + sigma_d * log_dur_raw);
  vector[P] a_start;
  vector[P] b_end;
  {
    int g = 1;
    for (s in 1:S) {
      int p0 = phase_start[s];
      mid[p0] = mid_first[s];
      // Gaps run from older to younger, so each subsequent midpoint is smaller.
      for (k in 1:(n_phase[s] - 1)) {
        mid[p0 + k] = mid[p0 + k - 1] - exp(log_gap[g]);
        g += 1;
      }
    }
  }
  a_start = mid + 0.5 * dur;
  b_end   = mid - 0.5 * dur;
}

model {
  // Priors
  // mid_first is uniform over the calendar range through its declared bounds.
  log_gap ~ normal(log(150), 1);      // lognormal separation between episodes
  log_dur_raw ~ std_normal();
  mu_d ~ normal(mu_d_prior_mean, 1);
  sigma_d ~ normal(0, 1);             // half-normal through the lower bound
  rho ~ beta(2, 18);
  for (m in 1:n_site_J2) pi2[m] ~ dirichlet(rep_vector(0.5, 2));   // alpha / J
  for (m in 1:n_site_J4) pi4[m] ~ dirichlet(rep_vector(0.25, 4));
  for (c in 1:2) q[c] ~ dirichlet(rep_vector(1.0, n_kappa));

  if (prior_only == 1) return;

  // Likelihood
  for (i in 1:N) {
    int s = site[i];
    int p0 = phase_start[s];
    int J = n_phase[s];
    int c = class_id[i];

    vector[J] lp_phase;
    for (j in 1:J) {
      int p = p0 + j - 1;
      vector[n_var[i]] lp_var;
      for (m in 1:n_var[i]) {
        int v = var_start[i] + m - 1;
        real lw;
        if (is_outlier[v] == 1) {
          lw = log(rho[c]);
        } else if (c == 1) {
          lw = log1m(rho[c]);                       // one kappa variant only
        } else {
          lw = log1m(rho[c]) + log(q[c - 1][kappa_idx[v]]);
        }
        real A = phase_average(a_start[p], b_end[p], v, L_flat, Phi_flat,
                               curve_pos, curve_len, curve_t0, phi_total, dt);
        lp_var[m] = lw + log(A + 1e-300);
      }
      real lpi;
      if (J == 1) lpi = 0;
      else if (J == 2) lpi = log(pi2[slot2[s]][j]);
      else lpi = log(pi4[slot4[s]][j]);
      lp_phase[j] = lpi + log_sum_exp(lp_var);
    }
    target += log_sum_exp(lp_phase);
  }
}

generated quantities {
  vector[N] log_lik;
  array[N] int z;                 // phase allocation, drawn from its conditional
  array[P] int<lower=0, upper=1> active = rep_array(0, P);
  array[P] real phase_start_calBP;
  array[P] real phase_end_calBP;
  array[P] real phase_duration;
  array[P] int n_allocated = rep_array(0, P);

  for (p in 1:P) {
    phase_start_calBP[p] = a_start[p];
    phase_end_calBP[p] = b_end[p];
    phase_duration[p] = dur[p];
  }

  for (i in 1:N) {
    int s = site[i];
    int p0 = phase_start[s];
    int J = n_phase[s];
    int c = class_id[i];

    vector[J] lp_phase;
    for (j in 1:J) {
      int p = p0 + j - 1;
      vector[n_var[i]] lp_var;
      for (m in 1:n_var[i]) {
        int v = var_start[i] + m - 1;
        real lw;
        if (is_outlier[v] == 1) lw = log(rho[c]);
        else if (c == 1) lw = log1m(rho[c]);
        else lw = log1m(rho[c]) + log(q[c - 1][kappa_idx[v]]);
        real A = phase_average(a_start[p], b_end[p], v, L_flat, Phi_flat,
                               curve_pos, curve_len, curve_t0, phi_total, dt);
        lp_var[m] = lw + log(A + 1e-300);
      }
      real lpi;
      if (J == 1) lpi = 0;
      else if (J == 2) lpi = log(pi2[slot2[s]][j]);
      else lpi = log(pi4[slot4[s]][j]);
      lp_phase[j] = lpi + log_sum_exp(lp_var);
    }
    log_lik[i] = log_sum_exp(lp_phase);

    // Recover the discrete allocation that the likelihood marginalised away.
    z[i] = p0 + categorical_logit_rng(lp_phase) - 1;
    n_allocated[z[i]] += 1;
  }

  // A phase is an occupation episode when at least one determination is
  // allocated to it. This is the definition used throughout: an episode exists
  // when there is dated evidence for it, with no weight threshold.
  for (p in 1:P) active[p] = n_allocated[p] > 0 ? 1 : 0;
}
```

- [ ] **Step 2: Write the failing test**

```r
# tests/testthat/test_susq_phases_single.R
source("../../R/susq_data.R")
source("../../R/susq_calibration.R")
source("../../R/susq_stan_data.R")
source("../../R/susq_fit.R")
library(posterior)

test_that("the model recovers a known single occupation phase", {
  # Simulate 12 short-lived determinations uniformly spread across a 60-year
  # window centred on 650 cal BP, with realistic errors, then check that the
  # posterior for the boundaries covers the truth.
  set.seed(20260813)
  cc <- intcal20_on_grid()
  truth <- list(a = 680, b = 620)
  t_true <- runif(12, truth$b, truth$a)
  c14 <- approx(cc$calBP, cc$c14, xout = t_true)$y
  err <- rep(30, 12)
  obs <- round(rnorm(12, c14, sqrt(err^2 + approx(cc$calBP, cc$err, t_true)$y^2)))

  sim <- data.frame(
    det_id = 1:12, site = "SIM", context = as.character(1:12),
    lab_no = paste0("SIM-", 1:12), c14_age = obs, c14_error = err,
    material_class = "short_lived", class_id = 1L, reference = "sim",
    pair_lab_no = NA_character_, pair_c14_age = NA_real_,
    pair_c14_error = NA_real_, pair_offset = NA_real_, preprocess = "none",
    stringsAsFactors = FALSE)

  v <- build_variants(sim, cc)
  sd <- build_stan_data(sim, v, data.frame(site = "SIM", n = 12L, J = 1L))
  fit <- fit_susq(sd, dates = sim, chains = 4, iter_warmup = 1000,
                  iter_sampling = 1000, refresh = 0)

  dg <- check_diagnostics(fit, strict = FALSE)
  expect_lt(dg$max_rhat, 1.01)
  expect_equal(dg$n_divergent, 0)

  s <- summarise_draws(fit$draws(c("phase_start_calBP", "phase_end_calBP")),
                       ~quantile(.x, c(0.025, 0.5, 0.975)))
  a_ci <- as.numeric(s[1, 2:4]); b_ci <- as.numeric(s[2, 2:4])
  expect_lt(a_ci[1], truth$a); expect_gt(a_ci[3], truth$a)
  expect_lt(b_ci[1], truth$b); expect_gt(b_ci[3], truth$b)
})

test_that("the prior predictive runs and produces plausible durations", {
  cc <- intcal20_on_grid()
  sim <- data.frame(det_id = 1L, site = "SIM", context = "1", lab_no = "SIM-1",
                    c14_age = 650, c14_error = 30, material_class = "short_lived",
                    class_id = 1L, reference = "sim", pair_lab_no = NA_character_,
                    pair_c14_age = NA_real_, pair_c14_error = NA_real_,
                    pair_offset = NA_real_, preprocess = "none",
                    stringsAsFactors = FALSE)
  v <- build_variants(sim, cc)
  sd <- build_stan_data(sim, v, data.frame(site = "SIM", n = 1L, J = 1L))
  sd$prior_only <- 1L
  fit <- fit_susq(sd, dates = sim, chains = 2, iter_warmup = 500,
                  iter_sampling = 500, refresh = 0)
  dur <- as.numeric(fit$draws("phase_duration", format = "matrix"))
  expect_gt(median(dur), 15); expect_lt(median(dur), 250)
  expect_lt(quantile(dur, 0.99), 5000)
})
```

- [ ] **Step 3: Write `R/susq_fit.R`**

```r
# R/susq_fit.R
#' Compile, initialise and fit the occupation models, and check the sampler.

suppressPackageStartupMessages({ library(cmdstanr); library(posterior) })

SUSQ_SEED <- 20260813L

#' Initialise each site's oldest phase midpoint near the calibrated dates that
#' belong to it. Without this, a phase can start far from all data, where every
#' phase-average is zero and the gradient carries no information.
init_from_data <- function(stan_data, dates, chains) {
  cc <- intcal20_on_grid()
  approx_cal <- approx(cc$c14, cc$calBP, xout = dates$c14_age, rule = 2)$y
  by_site <- tapply(approx_cal, stan_data$site, median)
  centre <- as.numeric(by_site[as.character(seq_len(stan_data$S))])
  centre[is.na(centre)] <- median(approx_cal)

  lapply(seq_len(chains), function(k) {
    list(
      mid_first  = pmin(pmax(centre + rnorm(stan_data$S, 0, 25),
                             stan_data$cal_min + 50), stan_data$cal_max - 50),
      log_gap    = rnorm(max(stan_data$P - stan_data$S, 0), log(150), 0.2),
      log_dur_raw = rnorm(stan_data$P, 0, 0.3),
      mu_d       = stan_data$mu_d_prior_mean,
      sigma_d    = 0.5,
      rho        = rep(0.1, 3)
    )
  })
}

fit_susq <- function(stan_data, dates,
                     model_file = "stan/susq_phases.stan",
                     chains = 4, iter_warmup = 2000, iter_sampling = 2000,
                     adapt_delta = 0.95, max_treedepth = 12,
                     seed = SUSQ_SEED, refresh = 200) {
  if (is.null(stan_data$prior_only)) stan_data$prior_only <- 0L
  payload <- stan_data[setdiff(names(stan_data), "site_names")]

  mod <- cmdstan_model(model_file, include_paths = "stan")
  mod$sample(
    data = payload,
    chains = chains, parallel_chains = min(chains, parallel::detectCores() - 1L),
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    adapt_delta = adapt_delta, max_treedepth = max_treedepth,
    seed = seed, refresh = refresh,
    init = init_from_data(stan_data, dates, chains)
  )
}

check_diagnostics <- function(fit, strict = TRUE) {
  pars <- c("mid", "dur", "mu_d", "sigma_d", "rho")
  dr <- fit$draws(pars)
  s <- posterior::summarise_draws(dr, "rhat", "ess_bulk", "ess_tail")
  sink(nullfile()); dg <- fit$diagnostic_summary(quiet = TRUE); sink()

  out <- data.frame(
    max_rhat = max(s$rhat, na.rm = TRUE),
    min_ess_bulk = min(s$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(s$ess_tail, na.rm = TRUE),
    n_divergent = sum(dg$num_divergent),
    n_max_treedepth = sum(dg$num_max_treedepth),
    ebfmi_min = min(dg$ebfmi)
  )
  if (strict) {
    if (out$max_rhat >= 1.01) stop("R-hat ", round(out$max_rhat, 4), " >= 1.01")
    if (out$min_ess_bulk < 400) stop("bulk ESS ", round(out$min_ess_bulk), " < 400")
    if (out$min_ess_tail < 400) stop("tail ESS ", round(out$min_ess_tail), " < 400")
    if (out$n_divergent > 0) stop(out$n_divergent, " divergent transitions")
    if (out$ebfmi_min < 0.3) stop("E-BFMI ", round(out$ebfmi_min, 3), " < 0.3")
  }
  out
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test_susq_phases_single.R")'`
Expected: PASS, 2 tests. Allow several minutes for compilation and sampling.

If the boundaries do not cover the truth, the likely cause is the sign convention on `mid`, `a_start` or `b_end`. Confirm that `a_start > b_end` and that both are in cal BP.

- [ ] **Step 5: Commit**

```bash
git add stan/susq_phases.stan R/susq_fit.R tests/testthat/test_susq_phases_single.R
git commit -m "Add multi-phase occupation model with single-phase recovery test"
```

---

### Task 6: Multi-phase recovery and the sequential model

**Files:**
- Create: `stan/susq_sequential.stan`
- Test: `tests/testthat/test_susq_multiphase.R`

**Interfaces:**
- Consumes: everything from Tasks 1 to 5.
- Produces: `stan/susq_sequential.stan`, the strictly non-overlapping model, taking the same data as `susq_phases.stan` plus nothing extra. It requires `n_phase` to be all 1 and uses a global ordered boundary vector. Sites must be supplied in ascending order of median calibrated date, which `build_stan_data(..., model = "sequential")` arranges.

- [ ] **Step 1: Write the sequential model**

```stan
// stan/susq_sequential.stan
/**
 * Strictly serial occupation. Every site has one phase and no two phases
 * overlap. Boundaries form a single ordered vector, so phase k ends before
 * phase k+1 begins.
 *
 * The order of sites is fixed by the data, ascending in median calibrated
 * date. That is the ordering most favourable to the serial hypothesis, chosen
 * using the data, so a win for this model is suggestive rather than decisive
 * while a loss is correspondingly strong. See spec section 5.
 */
functions {
#include susq_functions.stan
}

data {
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1> P;                                  // equals S here
  array[N] int<lower=1, upper=S> site;
  array[S] int<lower=1, upper=P> phase_start;
  array[S] int<lower=1> n_phase;
  array[N] int<lower=1, upper=3> class_id;

  int<lower=1> V;
  array[N] int<lower=1, upper=V> var_start;
  array[N] int<lower=1> n_var;
  array[V] int<lower=0> kappa_idx;
  array[V] int<lower=0, upper=1> is_outlier;

  int<lower=1> Ltot;
  vector[Ltot] L_flat;
  vector[Ltot] Phi_flat;
  array[V] int<lower=1> curve_pos;
  array[V] int<lower=1> curve_len;
  array[V] real curve_t0;
  array[V] real phi_total;
  real<lower=0> dt;

  int<lower=1> n_kappa;
  real cal_min;
  real cal_max;
  real mu_d_prior_mean;
  int<lower=0, upper=1> prior_only;
}

parameters {
  // 2S boundaries, increasing in cal BP: b_1 < a_1 < b_2 < a_2 < ...
  // so phase 1 is the youngest and phase S the oldest, and none overlap.
  // The calendar range is imposed by the declaration, which gives a proper
  // uniform prior over ordered boundaries inside it.
  ordered<lower=cal_min, upper=cal_max>[2 * S] bound;
  real mu_d;
  real<lower=0> sigma_d;
  array[2] simplex[n_kappa] q;
  vector<lower=0, upper=1>[3] rho;
}

transformed parameters {
  vector[S] b_end;
  vector[S] a_start;
  vector[S] dur;
  for (s in 1:S) {
    b_end[s] = bound[2 * s - 1];
    a_start[s] = bound[2 * s];
    dur[s] = a_start[s] - b_end[s];
  }
}

model {
  // bound is uniform over its declared ordered range. The map from
  // (b_end, a_start) to (b_end, dur) is linear with unit Jacobian, so the
  // duration prior below needs no adjustment.
  dur ~ lognormal(mu_d, sigma_d);
  mu_d ~ normal(mu_d_prior_mean, 1);
  sigma_d ~ normal(0, 1);
  rho ~ beta(2, 18);
  for (c in 1:2) q[c] ~ dirichlet(rep_vector(1.0, n_kappa));

  if (prior_only == 1) return;

  for (i in 1:N) {
    int s = site[i];
    int c = class_id[i];
    vector[n_var[i]] lp_var;
    for (m in 1:n_var[i]) {
      int v = var_start[i] + m - 1;
      real lw;
      if (is_outlier[v] == 1) lw = log(rho[c]);
      else if (c == 1) lw = log1m(rho[c]);
      else lw = log1m(rho[c]) + log(q[c - 1][kappa_idx[v]]);
      real A = phase_average(a_start[s], b_end[s], v, L_flat, Phi_flat,
                             curve_pos, curve_len, curve_t0, phi_total, dt);
      lp_var[m] = lw + log(A + 1e-300);
    }
    target += log_sum_exp(lp_var);
  }
}

generated quantities {
  vector[N] log_lik;
  array[S] real phase_start_calBP;
  array[S] real phase_end_calBP;
  array[S] real phase_duration;
  array[S] int<lower=0, upper=1> active = rep_array(1, S);

  for (s in 1:S) {
    phase_start_calBP[s] = a_start[s];
    phase_end_calBP[s] = b_end[s];
    phase_duration[s] = dur[s];
  }
  for (i in 1:N) {
    int s = site[i];
    int c = class_id[i];
    vector[n_var[i]] lp_var;
    for (m in 1:n_var[i]) {
      int v = var_start[i] + m - 1;
      real lw;
      if (is_outlier[v] == 1) lw = log(rho[c]);
      else if (c == 1) lw = log1m(rho[c]);
      else lw = log1m(rho[c]) + log(q[c - 1][kappa_idx[v]]);
      real A = phase_average(a_start[s], b_end[s], v, L_flat, Phi_flat,
                             curve_pos, curve_len, curve_t0, phi_total, dt);
      lp_var[m] = lw + log(A + 1e-300);
    }
    log_lik[i] = log_sum_exp(lp_var);
  }
}
```

- [ ] **Step 2: Extend `build_stan_data` for the sequential configuration**

Add `"sequential"` to the `model` argument in `R/susq_stan_data.R`. It behaves like `"indep"` except that sites are reordered by ascending median calibrated date before `site_id` is assigned. Insert this immediately after the `if (model == "single")` line:

```r
  if (model == "sequential") {
    J <- rep(1L, length(sites))
    cc <- intcal20_on_grid()
    approx_cal <- approx(cc$c14, cc$calBP, xout = dates$c14_age, rule = 2)$y
    med <- tapply(approx_cal, dates$site, median)
    sites <- names(sort(med))            # youngest first, matching the ordered vector
    site_id <- match(dates$site, sites)
  }
```

- [ ] **Step 3: Write the failing test**

```r
# tests/testthat/test_susq_multiphase.R
source("../../R/susq_data.R"); source("../../R/susq_calibration.R")
source("../../R/susq_stan_data.R"); source("../../R/susq_fit.R")
library(posterior)

simulate_site <- function(name, windows, n_each, err = 30, seed = 1) {
  set.seed(seed); cc <- intcal20_on_grid()
  rows <- list()
  for (w in seq_along(windows)) {
    tt <- runif(n_each[w], windows[[w]][2], windows[[w]][1])
    mu <- approx(cc$calBP, cc$c14, xout = tt)$y
    ce <- approx(cc$calBP, cc$err, xout = tt)$y
    rows[[w]] <- data.frame(site = name, c14_age = round(rnorm(length(tt), mu,
                              sqrt(err^2 + ce^2))), c14_error = err)
  }
  out <- do.call(rbind, rows)
  data.frame(det_id = seq_len(nrow(out)), site = out$site,
             context = as.character(seq_len(nrow(out))),
             lab_no = paste0(name, "-", seq_len(nrow(out))),
             c14_age = out$c14_age, c14_error = out$c14_error,
             material_class = "short_lived", class_id = 1L, reference = "sim",
             pair_lab_no = NA_character_, pair_c14_age = NA_real_,
             pair_c14_error = NA_real_, pair_offset = NA_real_,
             preprocess = "none", stringsAsFactors = FALSE)
}

test_that("two well-separated episodes at one site are both recovered", {
  sim <- simulate_site("SIM", list(c(1000, 950), c(500, 450)),
                       c(8, 8), seed = 20260813)
  cc <- intcal20_on_grid()
  v <- build_variants(sim, cc)
  sd <- build_stan_data(sim, v, data.frame(site = "SIM", n = 16L, J = 4L))
  fit <- fit_susq(sd, dates = sim, refresh = 0)

  expect_lt(check_diagnostics(fit, strict = FALSE)$max_rhat, 1.01)

  act <- fit$draws("active", format = "matrix")
  n_active <- rowSums(act)
  # The truth is two episodes. The posterior mode should be 2, and the model
  # must not routinely claim only one.
  expect_equal(as.integer(names(which.max(table(n_active)))), 2L)
  expect_lt(mean(n_active == 1), 0.15)

  # Both true windows should be covered by some active phase in most draws.
  st <- fit$draws("phase_start_calBP", format = "matrix")
  en <- fit$draws("phase_end_calBP", format = "matrix")
  covers <- function(t) mean(rowSums((st >= t) & (en <= t) & (act == 1)) > 0)
  expect_gt(covers(975), 0.8)
  expect_gt(covers(475), 0.8)
})

test_that("a single episode is not split into spurious extra episodes", {
  sim <- simulate_site("SIM", list(c(700, 650)), 12, seed = 7)
  cc <- intcal20_on_grid()
  v <- build_variants(sim, cc)
  sd <- build_stan_data(sim, v, data.frame(site = "SIM", n = 12L, J = 4L))
  fit <- fit_susq(sd, dates = sim, refresh = 0)
  n_active <- rowSums(fit$draws("active", format = "matrix"))
  # The overfitted-mixture prior should empty the unneeded components.
  expect_gt(mean(n_active <= 2), 0.8)
})

test_that("the sequential model compiles, runs and enforces non-overlap", {
  a <- simulate_site("A", list(c(900, 860)), 6, seed = 1)
  b <- simulate_site("B", list(c(500, 460)), 6, seed = 2); b$site <- "B"
  sim <- rbind(a, b); sim$det_id <- seq_len(nrow(sim))
  sim$lab_no <- paste0(sim$site, "-", sim$det_id)
  cc <- intcal20_on_grid()
  v <- build_variants(sim, cc)
  sd <- build_stan_data(sim, v, site_phase_counts(sim), model = "sequential")
  fit <- fit_susq(sd, dates = sim, model_file = "stan/susq_sequential.stan",
                  refresh = 0)
  st <- fit$draws("phase_start_calBP", format = "matrix")
  en <- fit$draws("phase_end_calBP", format = "matrix")
  # Phase 1 is the youngest: its start must be older than nothing else's end.
  expect_true(all(st[, 1] <= en[, 2] + 1e-8))
})
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test_susq_multiphase.R")'`
Expected: PASS, 3 tests.

If the two-episode test fails with the posterior collapsing to one phase, check the `log_gap` prior. A prior that is too tight around small gaps prevents the mixture from separating episodes.

- [ ] **Step 5: Commit**

```bash
git add stan/susq_sequential.stan R/susq_stan_data.R tests/testthat/test_susq_multiphase.R
git commit -m "Add sequential occupation model and multi-phase recovery tests"
```

---

### Task 7: Derived quantities in R

**Files:**
- Create: `R/susq_metrics.R`
- Test: `tests/testthat/test_susq_metrics.R`

**Interfaces:**
- Consumes: a fitted `CmdStanMCMC` and the `stan_data` used to fit it.
- Produces:
  - `occupancy_draws(fit, stan_data, report_dt = 10)` returns a list with `t` (reporting grid in cal BP) and `O` (a draws-by-sites-by-grid logical array), where a site is occupied at `t` if any active phase covers it.
  - `count_curve(occ, sites = NULL)` returns a data frame with `calBP`, `AD`, `median`, `lo50`, `hi50`, `lo95`, `hi95`, `mean`, from the posterior of the number of occupied sites. `sites` restricts the count to a named subset.
  - `pairwise_overlap(occ, site_names)` returns a data frame with `site_a`, `site_b`, `median_years`, `lo95`, `hi95`, `p_ge_25`, `p_any`.
  - `phase_summary(fit, stan_data)` returns one row per candidate phase with `site`, `phase`, `p_active`, and posterior quantiles of `start_calBP`, `end_calBP`, `duration`, conditional on the phase being active.

- [ ] **Step 1: Write the failing test**

```r
# tests/testthat/test_susq_metrics.R
source("../../R/susq_metrics.R")

# A hand-built fake draws object lets these be tested without sampling.
fake <- function() {
  # 2 sites, 1 phase each, 4 draws.
  # Site 1: 800-700 in every draw. Site 2: 750-650, but inactive in draw 4.
  list(
    stan_data = list(S = 2L, P = 2L, n_phase = c(1L, 1L),
                     phase_start = c(1L, 2L), site_names = c("A", "B")),
    start  = cbind(c(800, 800, 800, 800), c(750, 750, 750, 750)),
    end    = cbind(c(700, 700, 700, 700), c(650, 650, 650, 650)),
    active = cbind(c(1, 1, 1, 1), c(1, 1, 1, 0))
  )
}

test_that("occupancy marks a site occupied only where an active phase covers t", {
  f <- fake()
  occ <- occupancy_from_arrays(f$start, f$end, f$active, f$stan_data,
                               report_dt = 10)
  i760 <- which(occ$t == 760); i690 <- which(occ$t == 690)
  expect_true(all(occ$O[, 1, i760]))               # A covers 760 always
  expect_false(any(occ$O[, 1, which(occ$t == 660)]))  # A ends at 700
  expect_equal(sum(occ$O[, 2, i690]), 3)           # B active in 3 of 4 draws
})

test_that("count_curve counts sites, not phases, and reports credible bands", {
  f <- fake()
  occ <- occupancy_from_arrays(f$start, f$end, f$active, f$stan_data, 10)
  cc <- count_curve(occ)
  expect_true(all(c("calBP","AD","median","lo95","hi95") %in% names(cc)))
  expect_equal(cc$AD, 1950 - cc$calBP)
  at740 <- cc[cc$calBP == 740, ]
  expect_equal(at740$median, 2)                    # both sites cover 740
  at660 <- cc[cc$calBP == 660, ]
  expect_equal(at660$median, 1)                    # only B
  expect_true(all(cc$lo95 <= cc$median & cc$median <= cc$hi95))
})

test_that("count_curve restricted to a subset counts only those sites", {
  f <- fake()
  occ <- occupancy_from_arrays(f$start, f$end, f$active, f$stan_data, 10)
  cc <- count_curve(occ, sites = "A")
  expect_true(all(cc$hi95 <= 1))
})

test_that("pairwise_overlap reports years, not a saturating indicator", {
  f <- fake()
  occ <- occupancy_from_arrays(f$start, f$end, f$active, f$stan_data, 10)
  po <- pairwise_overlap(occ, f$stan_data$site_names)
  expect_equal(nrow(po), 1L)
  # A is 800-700, B is 750-650, so the overlap is 700-750, about 50 years,
  # in the 3 draws where B is active.
  expect_gt(po$median_years, 40); expect_lt(po$median_years, 60)
  expect_equal(po$p_any, 0.75)
  expect_equal(po$p_ge_25, 0.75)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test_susq_metrics.R")'`
Expected: FAIL, `could not find function "occupancy_from_arrays"`.

- [ ] **Step 3: Write the implementation**

```r
# R/susq_metrics.R
#' Derived quantities computed in R from the posterior draws.
#' Keeping these out of Stan keeps the fit output small and lets the same draws
#' answer questions posed after the fact, such as restricting the count to
#' well-dated sites. See spec section 4.

suppressPackageStartupMessages({ library(posterior) })

#' Site-by-calendar-year occupancy for every draw.
#' A site is occupied at t when at least one of its active phases covers t.
occupancy_from_arrays <- function(start, end, active, stan_data, report_dt = 10) {
  t_grid <- seq(0, 2000, by = report_dt)
  D <- nrow(start); S <- stan_data$S; G <- length(t_grid)
  O <- array(FALSE, dim = c(D, S, G))
  for (s in seq_len(S)) {
    p_idx <- stan_data$phase_start[s] + seq_len(stan_data$n_phase[s]) - 1L
    for (p in p_idx) {
      cover <- outer(start[, p], t_grid, ">=") & outer(end[, p], t_grid, "<=")
      cover <- cover & (active[, p] == 1)
      O[, s, ] <- O[, s, ] | cover
    }
  }
  list(t = t_grid, O = O, sites = stan_data$site_names, dt = report_dt)
}

occupancy_draws <- function(fit, stan_data, report_dt = 10) {
  occupancy_from_arrays(
    fit$draws("phase_start_calBP", format = "matrix"),
    fit$draws("phase_end_calBP", format = "matrix"),
    fit$draws("active", format = "matrix"),
    stan_data, report_dt)
}

#' Posterior number of simultaneously occupied sites at each calendar year.
count_curve <- function(occ, sites = NULL) {
  keep <- if (is.null(sites)) seq_along(occ$sites) else which(occ$sites %in% sites)
  counts <- apply(occ$O[, keep, , drop = FALSE], c(1, 3), sum)
  q <- function(p) apply(counts, 2, quantile, probs = p, names = FALSE)
  data.frame(
    calBP = occ$t, AD = 1950 - occ$t,
    mean = colMeans(counts),
    median = q(0.5),
    lo50 = q(0.25), hi50 = q(0.75),
    lo95 = q(0.025), hi95 = q(0.975)
  )
}

#' Years of joint occupation for every pair of sites.
pairwise_overlap <- function(occ, site_names) {
  S <- length(site_names); rows <- list(); k <- 0L
  for (i in seq_len(S - 1)) for (j in (i + 1):S) {
    years <- rowSums(occ$O[, i, ] & occ$O[, j, ]) * occ$dt
    k <- k + 1L
    rows[[k]] <- data.frame(
      site_a = site_names[i], site_b = site_names[j],
      median_years = median(years),
      lo95 = unname(quantile(years, 0.025)),
      hi95 = unname(quantile(years, 0.975)),
      p_ge_25 = mean(years >= 25),
      p_any = mean(years > 0),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

#' Per-phase summary, with boundary and duration quantiles conditional on the
#' phase being active. An inactive phase has no archaeological meaning, so
#' summarising it unconditionally would mix real episodes with empty ones.
phase_summary <- function(fit, stan_data) {
  st <- fit$draws("phase_start_calBP", format = "matrix")
  en <- fit$draws("phase_end_calBP", format = "matrix")
  du <- fit$draws("phase_duration", format = "matrix")
  ac <- fit$draws("active", format = "matrix")

  rows <- list(); k <- 0L
  for (s in seq_len(stan_data$S)) {
    for (jj in seq_len(stan_data$n_phase[s])) {
      p <- stan_data$phase_start[s] + jj - 1L
      sel <- ac[, p] == 1
      k <- k + 1L
      qs <- function(v) if (sum(sel) < 20) rep(NA_real_, 3) else
        unname(quantile(v[sel], c(0.025, 0.5, 0.975)))
      rows[[k]] <- data.frame(
        site = stan_data$site_names[s], phase = jj,
        p_active = mean(sel),
        start_lo = qs(st[, p])[1], start_med = qs(st[, p])[2], start_hi = qs(st[, p])[3],
        end_lo = qs(en[, p])[1], end_med = qs(en[, p])[2], end_hi = qs(en[, p])[3],
        dur_lo = qs(du[, p])[1], dur_med = qs(du[, p])[2], dur_hi = qs(du[, p])[3],
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test_susq_metrics.R")'`
Expected: PASS, 4 tests.

- [ ] **Step 5: Commit**

```bash
git add R/susq_metrics.R tests/testthat/test_susq_metrics.R
git commit -m "Add occupancy, contemporaneous count, pairwise overlap and phase summaries"
```

---

### Task 8: Simulation recovery

**Files:**
- Create: `R/susq_simulate.R`
- Create: `scripts/03_simulation_recovery.R`
- Test: `tests/testthat/test_susq_simulate.R`

**Interfaces:**
- Consumes: Tasks 1 to 7.
- Produces:
  - `simulate_assemblage(dates, truth_seed)` returns a list with `dates` (same columns as `prepare_susq_data()$dates`, same site sizes, errors and material classes as the real data) and `truth` (data frame of `site`, `phase`, `start_calBP`, `end_calBP`).
  - `recovery_run(rep_id)` fits one replicate and returns coverage indicators.
  - `scripts/03_simulation_recovery.R` runs 30 replicates and writes `output/susquehanna/simulation_recovery.rds` plus a coverage table.

This task is what determines which sites can carry an interpretation. The expectation stated in the spec is that multi-phase structure will not be recoverable at three determinations per site. The point is to measure that, not to assume it.

- [ ] **Step 1: Write the failing test**

```r
# tests/testthat/test_susq_simulate.R
source("../../R/susq_data.R"); source("../../R/susq_calibration.R")
source("../../R/susq_simulate.R")

test_that("simulated assemblages match the real one in shape", {
  real <- prepare_susq_data()$dates
  sim <- simulate_assemblage(real, truth_seed = 1)
  expect_equal(nrow(sim$dates), nrow(real))
  expect_equal(sort(table(sim$dates$site)), sort(table(real$site)))
  expect_equal(sort(table(sim$dates$material_class)),
               sort(table(real$material_class)))
  expect_equal(sim$dates$c14_error, real$c14_error)
  expect_true(all(sim$truth$start_calBP > sim$truth$end_calBP))
})

test_that("wood determinations are simulated with inbuilt age and short-lived are not", {
  real <- prepare_susq_data()$dates
  sim <- simulate_assemblage(real, truth_seed = 2)
  expect_true(!is.null(sim$inbuilt))
  expect_true(all(sim$inbuilt[sim$dates$material_class == "short_lived"] == 0))
  expect_true(mean(sim$inbuilt[sim$dates$material_class == "wood"]) > 0)
})

test_that("simulated sites with more phases get more phases in the truth", {
  real <- prepare_susq_data()$dates
  sim <- simulate_assemblage(real, truth_seed = 3)
  pc <- site_phase_counts(real)
  n_truth <- table(sim$truth$site)
  big <- pc$site[pc$J == 4L]
  expect_true(all(n_truth[big] >= 1))
  expect_true(all(n_truth <= pc$J[match(names(n_truth), pc$site)]))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test_susq_simulate.R")'`
Expected: FAIL, `could not find function "simulate_assemblage"`.

- [ ] **Step 3: Write the implementation**

```r
# R/susq_simulate.R
#' Generate synthetic assemblages that match the real one in every respect
#' except that the occupation history is known. See spec section 6.2.

simulate_assemblage <- function(dates, truth_seed) {
  set.seed(truth_seed)
  cc <- intcal20_on_grid()
  pc <- site_phase_counts(dates)

  truth <- list(); k <- 0L
  phase_of <- integer(nrow(dates))
  t_event <- numeric(nrow(dates))

  for (s in seq_len(nrow(pc))) {
    site <- pc$site[s]
    rows <- which(dates$site == site)
    # Draw the true number of episodes: one for small sites, one or two for
    # mid-sized, one to three for well-dated sites. Never more than J.
    n_ep <- min(pc$J[s], sample(seq_len(min(3L, pc$J[s])), 1))
    mids <- sort(runif(n_ep, 250, 1500), decreasing = TRUE)
    durs <- rlnorm(n_ep, log(60), 0.6)
    for (e in seq_len(n_ep)) {
      k <- k + 1L
      truth[[k]] <- data.frame(site = site, phase = e,
                               start_calBP = mids[e] + durs[e] / 2,
                               end_calBP = mids[e] - durs[e] / 2,
                               stringsAsFactors = FALSE)
    }
    assign_ep <- sample(seq_len(n_ep), length(rows), replace = TRUE)
    phase_of[rows] <- assign_ep
    for (r in seq_along(rows)) {
      e <- assign_ep[r]
      t_event[rows[r]] <- runif(1, mids[e] - durs[e] / 2, mids[e] + durs[e] / 2)
    }
  }
  truth <- do.call(rbind, truth)

  # Inbuilt age: exponential with mean 40 for wood, 20 for indeterminate,
  # 0 for short-lived. The sample's growth date is older than the event.
  mean_inbuilt <- c(short_lived = 0, wood = 40, indeterminate = 20)
  inbuilt <- vapply(dates$material_class, function(cl) {
    m <- mean_inbuilt[[cl]]
    if (m == 0) 0 else rexp(1, rate = 1 / m)
  }, 0)
  t_sample <- t_event + inbuilt

  mu <- approx(cc$calBP, cc$c14, xout = t_sample, rule = 2)$y
  ce <- approx(cc$calBP, cc$err, xout = t_sample, rule = 2)$y
  sim <- dates
  sim$c14_age <- round(rnorm(nrow(dates), mu, sqrt(dates$c14_error^2 + ce^2)))
  # Simulated data has no replicate or wiggle-match structure.
  sim$pair_lab_no <- NA_character_; sim$pair_c14_age <- NA_real_
  sim$pair_c14_error <- NA_real_; sim$pair_offset <- NA_real_
  sim$preprocess <- "none"

  list(dates = sim, truth = truth, inbuilt = inbuilt,
       t_event = t_event, phase_of = phase_of)
}

#' Fit one simulated assemblage and report whether the truth was covered.
recovery_run <- function(rep_id, real_dates) {
  sim <- simulate_assemblage(real_dates, truth_seed = 1000L + rep_id)
  cc <- intcal20_on_grid()
  v <- build_variants(sim$dates, cc)
  sd <- build_stan_data(sim$dates, v, site_phase_counts(sim$dates))
  fit <- fit_susq(sd, dates = sim$dates, refresh = 0)

  dg <- check_diagnostics(fit, strict = FALSE)
  occ <- occupancy_draws(fit, sd)
  cnt <- count_curve(occ)

  # True count at each reporting year.
  true_count <- vapply(occ$t, function(t)
    sum(vapply(split(sim$truth, sim$truth$site), function(x)
      any(x$start_calBP >= t & x$end_calBP <= t), TRUE)), 0)

  ps <- phase_summary(fit, sd)
  n_truth <- table(sim$truth$site)
  n_est <- tapply(ps$p_active, ps$site, sum)

  list(rep_id = rep_id, diagnostics = dg,
       count_coverage = mean(true_count >= cnt$lo95 & true_count <= cnt$hi95),
       count_bias = mean(cnt$median - true_count),
       n_phase_truth = as.integer(n_truth[names(n_est)]),
       n_phase_est = as.numeric(n_est),
       site = names(n_est),
       n_dates = as.integer(table(sim$dates$site)[names(n_est)]))
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test_susq_simulate.R")'`
Expected: PASS, 3 tests.

- [ ] **Step 5: Write the recovery script**

```r
# scripts/03_simulation_recovery.R
#' Run 30 simulation replicates and report coverage. See spec section 6.2.
for (f in c("R/susq_data.R", "R/susq_calibration.R", "R/susq_stan_data.R",
            "R/susq_fit.R", "R/susq_metrics.R", "R/susq_simulate.R")) source(f)
dir.create("output/susquehanna", showWarnings = FALSE, recursive = TRUE)

real <- prepare_susq_data()$dates
res <- lapply(1:30, function(i) { cat("replicate", i, "\n"); recovery_run(i, real) })
saveRDS(res, "output/susquehanna/simulation_recovery.rds")

cov <- vapply(res, `[[`, 0, "count_coverage")
bias <- vapply(res, `[[`, 0, "count_bias")
cat(sprintf("\nCount curve 95%% coverage: %.3f (target 0.95)\n", mean(cov)))
cat(sprintf("Count curve median bias: %+.2f sites\n", mean(bias)))

by_size <- do.call(rbind, lapply(res, function(r)
  data.frame(site = r$site, n = r$n_dates,
             truth = r$n_phase_truth, est = r$n_phase_est)))
by_size$bucket <- cut(by_size$n, c(0, 2, 5, 100),
                      labels = c("1-2 dates", "3-5 dates", "6+ dates"))
agg <- aggregate(cbind(truth, est) ~ bucket, by_size, mean)
agg$error <- agg$est - agg$truth
cat("\nEpisode-count recovery by site size:\n"); print(agg)
write.csv(agg, "output/susquehanna/recovery_by_site_size.csv", row.names = FALSE)
```

- [ ] **Step 6: Run it**

Run: `Rscript scripts/03_simulation_recovery.R 2>&1 | tail -30`
Expected: coverage near 0.95. If coverage is far below 0.95, the model is overconfident and the cause must be found before proceeding. Record the by-size table; it determines which sites carry interpretation.

- [ ] **Step 7: Commit**

```bash
git add R/susq_simulate.R scripts/03_simulation_recovery.R tests/testthat/test_susq_simulate.R output/susquehanna/recovery_by_site_size.csv
git commit -m "Add simulation recovery framework and coverage report"
```

---

### Task 9: Fit the real data and compare models

**Files:**
- Create: `scripts/01_prepare_data.R`
- Create: `scripts/02_prior_predictive.R`
- Create: `scripts/04_fit_models.R`
- Create: `scripts/05_compare_models.R`

**Interfaces:**
- Produces `output/susquehanna/fits/<name>.rds` for `multi`, `indep`, `single`, `sequential`, `multi_dur25`, `multi_dur150`, and `output/susquehanna/model_comparison.csv`.

- [ ] **Step 1: Write the data preparation script**

```r
# scripts/01_prepare_data.R
for (f in c("R/susq_data.R", "R/susq_calibration.R", "R/susq_stan_data.R")) source(f)
dir.create("output/susquehanna", showWarnings = FALSE, recursive = TRUE)

d <- prepare_susq_data()
cc <- intcal20_on_grid()
v <- build_variants(d$dates, cc)
pc <- site_phase_counts(d$dates)

saveRDS(list(dates = d$dates, quality = d$quality, curve = cc,
             variants = v, phase_counts = pc),
        "output/susquehanna/prepared.rds")

write.csv(d$dates, "output/susquehanna/determinations.csv", row.names = FALSE)
write.csv(d$quality, "output/susquehanna/data_quality.csv", row.names = FALSE)
write.csv(pc, "output/susquehanna/site_phase_counts.csv", row.names = FALSE)

cat(sprintf("%d determinations, %d sites, %d candidate phases\n",
            nrow(d$dates), nrow(pc), sum(pc$J)))
cat("\nQuality flags:\n"); print(d$quality)
```

- [ ] **Step 2: Write the prior predictive script**

```r
# scripts/02_prior_predictive.R
for (f in c("R/susq_data.R","R/susq_calibration.R","R/susq_stan_data.R",
            "R/susq_fit.R","R/susq_metrics.R")) source(f)
p <- readRDS("output/susquehanna/prepared.rds")

sd <- build_stan_data(p$dates, p$variants, p$phase_counts)
sd$prior_only <- 1L
fit <- fit_susq(sd, dates = p$dates, iter_warmup = 1000, iter_sampling = 1000)
saveRDS(fit$draws(), "output/susquehanna/prior_predictive_draws.rds")

dur <- as.numeric(fit$draws("phase_duration", format = "matrix"))
occ <- occupancy_draws(fit, sd); cnt <- count_curve(occ)
cat(sprintf("Prior duration: median %.0f, 95%% [%.0f, %.0f] years\n",
            median(dur), quantile(dur, 0.025), quantile(dur, 0.975)))
cat(sprintf("Prior count: median %.1f, 95%% up to %.0f sites\n",
            median(cnt$median), max(cnt$hi95)))
cat("\nCheck: durations should not routinely exceed a few centuries,\n",
    "and the prior count should not pin at 34 across the whole range.\n")
```

- [ ] **Step 3: Write the fitting script**

```r
# scripts/04_fit_models.R
for (f in c("R/susq_data.R","R/susq_calibration.R","R/susq_stan_data.R",
            "R/susq_fit.R","R/susq_metrics.R")) source(f)
p <- readRDS("output/susquehanna/prepared.rds")
dir.create("output/susquehanna/fits", showWarnings = FALSE, recursive = TRUE)

configs <- list(
  multi        = list(model = "multi",      mu_d = log(60),  file = "stan/susq_phases.stan"),
  indep        = list(model = "indep",      mu_d = log(60),  file = "stan/susq_phases.stan"),
  single       = list(model = "single",     mu_d = log(60),  file = "stan/susq_phases.stan"),
  sequential   = list(model = "sequential", mu_d = log(60),  file = "stan/susq_sequential.stan"),
  multi_dur25  = list(model = "multi",      mu_d = log(25),  file = "stan/susq_phases.stan"),
  multi_dur150 = list(model = "multi",      mu_d = log(150), file = "stan/susq_phases.stan")
)

diagnostics <- list()
for (nm in names(configs)) {
  cfg <- configs[[nm]]
  cat("\n=== fitting", nm, "===\n")
  sd <- build_stan_data(p$dates, p$variants, p$phase_counts,
                        mu_d_prior = cfg$mu_d, model = cfg$model)
  fit <- fit_susq(sd, dates = p$dates, model_file = cfg$file)
  dg <- check_diagnostics(fit, strict = FALSE)
  print(dg)
  diagnostics[[nm]] <- cbind(model = nm, dg)
  fit$save_object(sprintf("output/susquehanna/fits/%s.rds", nm))
  saveRDS(sd, sprintf("output/susquehanna/fits/%s_standata.rds", nm))
}
dg <- do.call(rbind, diagnostics)
write.csv(dg, "output/susquehanna/diagnostics.csv", row.names = FALSE)

bad <- dg[dg$max_rhat >= 1.01 | dg$n_divergent > 0 |
          dg$min_ess_bulk < 400 | dg$ebfmi_min < 0.3, ]
if (nrow(bad)) {
  cat("\nDIAGNOSTIC FAILURES. Do not interpret these fits:\n"); print(bad)
} else cat("\nAll fits pass the spec thresholds.\n")
```

- [ ] **Step 4: Write the comparison script**

```r
# scripts/05_compare_models.R
library(loo)
fits <- c("multi", "indep", "single", "sequential")
ll <- lapply(fits, function(nm)
  readRDS(sprintf("output/susquehanna/fits/%s.rds", nm))$draws("log_lik", format = "matrix"))
names(ll) <- fits

loos <- lapply(ll, function(m) loo(m, r_eff = relative_eff(exp(m), chain_id = rep(1:4, each = nrow(m) / 4))))
k_bad <- vapply(loos, function(l) mean(l$diagnostics$pareto_k > 0.7), 0)
print(data.frame(model = fits, frac_pareto_k_gt_0.7 = round(k_bad, 3)))

if (any(k_bad > 0.05)) {
  cat("\nPSIS-LOO is unreliable here (spec section 5). Run exact 10-fold CV\n",
      "with scripts/05b_kfold.R and report that instead.\n")
}
cmp <- loo_compare(loos)
print(cmp)
write.csv(as.data.frame(cmp), "output/susquehanna/model_comparison.csv")
```

- [ ] **Step 5: Run everything and record the outcome**

```bash
export PATH="/opt/homebrew/opt/gdal/bin:/opt/homebrew/opt/geos/bin:/opt/homebrew/opt/proj/bin:$PATH"
Rscript scripts/01_prepare_data.R
Rscript scripts/02_prior_predictive.R
Rscript scripts/04_fit_models.R
Rscript scripts/05_compare_models.R
```

Expected: all six fits pass R-hat < 1.01, ESS > 400, zero divergences. If any fail, fix before proceeding. Do not report results from a fit that failed diagnostics.

- [ ] **Step 6: Commit**

```bash
git add scripts/01_prepare_data.R scripts/02_prior_predictive.R scripts/04_fit_models.R scripts/05_compare_models.R output/susquehanna/diagnostics.csv output/susquehanna/model_comparison.csv output/susquehanna/data_quality.csv output/susquehanna/site_phase_counts.csv
git commit -m "Fit all six model configurations and compare by LOO"
```

---

### Task 10: Figures

**Files:**
- Create: `R/susq_theme.R`
- Create: `R/susq_plots.R`
- Create: `scripts/07_figures.R`
- Test: `tests/testthat/test_susq_theme.R`

**Interfaces:**
- Produces eight figures in `output/susquehanna/figures/`, named `fig1_count_through_time.png` through `fig8_diagnostics.png`.
- `theme_susq()` returns a ggplot theme. `OKABE_ITO` is the categorical palette. `susq_save(plot, file, width = 7, height = 4.5)` writes PNG at 300 dpi through `ragg`.

- [ ] **Step 1: Write the theme and its test**

```r
# R/susq_theme.R
#' Shared figure style. No titles; descriptive information belongs in captions.
suppressPackageStartupMessages({ library(ggplot2); library(systemfonts) })

OKABE_ITO <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

susq_font <- function() {
  fonts <- systemfonts::system_fonts()$family
  if ("Arial" %in% fonts) "Arial" else "sans"
}

theme_susq <- function(base_size = 10) {
  theme_minimal(base_size = base_size, base_family = susq_font()) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
      axis.line = element_line(linewidth = 0.35, colour = "grey20"),
      axis.ticks = element_line(linewidth = 0.35, colour = "grey20"),
      plot.title = element_blank(), plot.subtitle = element_blank(),
      legend.position = "bottom", legend.key.height = unit(0.8, "lines"),
      strip.text = element_text(face = "plain", hjust = 0)
    )
}

susq_save <- function(plot, file, width = 7, height = 4.5) {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  ggsave(file, plot, device = ragg::agg_png, width = width, height = height,
         units = "in", dpi = 300, bg = "white")
  invisible(file)
}

#' Calibration reversals measured from IntCal20, marked on time-axis figures so
#' that a change in interval width can be read against curve structure.
CURVE_REVERSALS <- data.frame(
  start_calBP = c(411, 886, 613, 1233, 1463, 1041),
  end_calBP   = c(344, 837, 577, 1199, 1429, 1014)
)
```

```r
# tests/testthat/test_susq_theme.R
source("../../R/susq_theme.R")
test_that("theme carries no title and the palette is Okabe-Ito", {
  th <- theme_susq()
  expect_s3_class(th, "theme")
  expect_s3_class(th$plot.title, "element_blank")
  expect_equal(length(OKABE_ITO), 8L)
  expect_equal(OKABE_ITO[2], "#E69F00")
})
test_that("susq_save writes a 300 dpi PNG of the requested width", {
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) +
    ggplot2::geom_point() + theme_susq()
  f <- tempfile(fileext = ".png")
  susq_save(p, f, width = 7, height = 4.5)
  expect_true(file.exists(f))
  info <- png::readPNG(f, info = TRUE)
  expect_equal(dim(info)[2], 7 * 300, tolerance = 2)
})
```

- [ ] **Step 2: Run the theme test**

Run: `Rscript -e 'install.packages("png", quiet=TRUE); testthat::test_file("tests/testthat/test_susq_theme.R")'`
Expected: PASS, 2 tests.

- [ ] **Step 3: Write the plotting functions**

```r
# R/susq_plots.R
#' Figures for the contemporaneity analysis. See spec section 7.
source("R/susq_theme.R")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

#' Figure 1. Posterior number of simultaneously occupied sites through time,
#' with the IntCal20 reversals marked.
plot_count_through_time <- function(cnt, cnt_subset = NULL) {
  rev_ad <- data.frame(xmin = 1950 - CURVE_REVERSALS$start_calBP,
                       xmax = 1950 - CURVE_REVERSALS$end_calBP)
  p <- ggplot(cnt, aes(AD)) +
    geom_rect(data = rev_ad, inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "grey90") +
    geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = OKABE_ITO[6], alpha = 0.20) +
    geom_ribbon(aes(ymin = lo50, ymax = hi50), fill = OKABE_ITO[6], alpha = 0.35) +
    geom_line(aes(y = median), colour = OKABE_ITO[6], linewidth = 0.7)
  if (!is.null(cnt_subset)) {
    p <- p + geom_line(data = cnt_subset, aes(y = median),
                       colour = OKABE_ITO[7], linewidth = 0.7, linetype = "22")
  }
  p + scale_x_continuous(breaks = seq(200, 1800, 200)) +
    labs(x = "Year AD", y = "Simultaneously occupied sites") +
    theme_susq()
}

#' Figure 2. Occupation intervals, one row per site, ordered by median start.
plot_occupation_intervals <- function(ps, dates, curve) {
  act <- ps[ps$p_active > 0.1 & !is.na(ps$start_med), ]
  act$label <- paste0(act$site, ifelse(act$phase > 1, paste0(" (", act$phase, ")"), ""))
  act <- act[order(act$start_med), ]
  act$label <- factor(act$label, levels = act$label)
  ggplot(act, aes(y = label)) +
    geom_segment(aes(x = 1950 - start_lo, xend = 1950 - end_hi,
                     yend = label), colour = "grey75", linewidth = 1.2) +
    geom_segment(aes(x = 1950 - start_med, xend = 1950 - end_med,
                     yend = label, alpha = p_active),
                 colour = OKABE_ITO[6], linewidth = 3) +
    scale_alpha_continuous(range = c(0.25, 1), name = "P(episode)") +
    labs(x = "Year AD", y = NULL) + theme_susq()
}

#' Figure 3. Pairwise overlap in years, and the probability of at least 25
#' years of joint occupation.
plot_overlap_matrix <- function(po, order_sites) {
  po2 <- rbind(po, transform(po, site_a = po$site_b, site_b = po$site_a))
  po2$site_a <- factor(po2$site_a, levels = order_sites)
  po2$site_b <- factor(po2$site_b, levels = rev(order_sites))
  ggplot(po2, aes(site_a, site_b, fill = median_years)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    scale_fill_viridis_c(name = "Median overlap (years)", na.value = "grey95") +
    labs(x = NULL, y = NULL) + theme_susq() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6), panel.grid = element_blank())
}

#' Figure 4. Duration posteriors against the hierarchical prior, so it is
#' visible which sites are data-driven and which are prior-driven.
plot_duration_posteriors <- function(fit, ps, mu_d_prior = log(60)) {
  du <- fit$draws("phase_duration", format = "matrix")
  ac <- fit$draws("active", format = "matrix")
  keep <- which(colMeans(ac) > 0.1)
  long <- do.call(rbind, lapply(keep, function(p) {
    sel <- ac[, p] == 1
    if (sum(sel) < 20) return(NULL)
    data.frame(phase = p, duration = du[sel, p],
               label = paste0(ps$site[p], if (ps$phase[p] > 1) paste0(" (", ps$phase[p], ")") else ""))
  }))
  ord <- tapply(long$duration, long$label, median)
  long$label <- factor(long$label, levels = names(sort(ord)))
  prior <- data.frame(duration = rlnorm(20000, mu_d_prior, 0.6))
  ggplot(long, aes(duration, label)) +
    ggridges::geom_density_ridges(scale = 1.6, fill = OKABE_ITO[4],
                                  colour = NA, alpha = 0.75, rel_min_height = 0.01) +
    geom_vline(xintercept = quantile(prior$duration, c(0.025, 0.5, 0.975)),
               linetype = c("22", "solid", "22"), colour = "grey40", linewidth = 0.3) +
    scale_x_continuous(trans = "log10", limits = c(3, 3000)) +
    labs(x = "Occupation duration (years)", y = NULL) + theme_susq()
}

#' Figure 5. The underlying evidence: calibrated determinations by site and
#' material class.
plot_calibrated_dates <- function(dates, curve) {
  rows <- lapply(seq_len(nrow(dates)), function(i) {
    L <- calib_likelihood(dates$c14_age[i], dates$c14_error[i], curve)
    L <- L / max(L)
    keep <- L > 0.01
    data.frame(det_id = dates$det_id[i], site = dates$site[i],
               material_class = dates$material_class[i],
               calBP = curve$calBP[keep], dens = L[keep])
  })
  df <- do.call(rbind, rows)
  df$AD <- 1950 - df$calBP
  ggplot(df, aes(AD, factor(det_id), height = dens, fill = material_class)) +
    ggridges::geom_ridgeline(scale = 2.5, colour = NA, alpha = 0.8) +
    facet_grid(rows = vars(site), scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_manual(values = OKABE_ITO[c(4, 7, 3)], name = NULL) +
    labs(x = "Year AD", y = NULL) + theme_susq() +
    theme(axis.text.y = element_blank(),
          strip.text.y.left = element_text(angle = 0, size = 5, hjust = 1))
}

#' Figure 6. LOO comparison.
plot_model_comparison <- function(cmp) {
  df <- as.data.frame(cmp); df$model <- rownames(df)
  df$model <- factor(df$model, levels = rev(df$model))
  ggplot(df, aes(elpd_diff, model)) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = elpd_diff - 2 * se_diff,
                       xmax = elpd_diff + 2 * se_diff), height = 0.15) +
    geom_point(size = 2.2, colour = OKABE_ITO[6]) +
    labs(x = "ELPD difference from the best model", y = NULL) + theme_susq()
}

#' Figure 7. Inbuilt-age posteriors against their priors. If a posterior tracks
#' its prior, the correction is an assumption and not an inference.
plot_inbuilt_age <- function(fit) {
  q <- fit$draws("q", format = "matrix")
  kappas <- c(0, 15, 30, 60, 120, 240)
  cls <- c("wood", "indeterminate")
  df <- do.call(rbind, lapply(1:2, function(c) do.call(rbind, lapply(1:6, function(k) {
    col <- sprintf("q[%d,%d]", c, k)
    data.frame(class = cls[c], kappa = kappas[k],
               med = median(q[, col]),
               lo = quantile(q[, col], 0.025), hi = quantile(q[, col], 0.975))
  }))))
  ggplot(df, aes(factor(kappa), med, colour = class)) +
    geom_hline(yintercept = 1 / 6, linetype = "22", colour = "grey50", linewidth = 0.3) +
    geom_pointrange(aes(ymin = lo, ymax = hi),
                    position = position_dodge(width = 0.4), size = 0.3) +
    scale_colour_manual(values = OKABE_ITO[c(7, 3)], name = NULL) +
    labs(x = "Inbuilt age (years)",
         y = "Posterior weight") + theme_susq()
}
```

- [ ] **Step 4: Write the figure script**

```r
# scripts/07_figures.R
for (f in c("R/susq_data.R","R/susq_calibration.R","R/susq_stan_data.R",
            "R/susq_fit.R","R/susq_metrics.R","R/susq_theme.R","R/susq_plots.R")) source(f)
library(loo)
p  <- readRDS("output/susquehanna/prepared.rds")
fit <- readRDS("output/susquehanna/fits/multi.rds")
sd  <- readRDS("output/susquehanna/fits/multi_standata.rds")
figs <- "output/susquehanna/figures"

occ <- occupancy_draws(fit, sd)
cnt <- count_curve(occ)
well_dated <- p$phase_counts$site[p$phase_counts$n >= 3]
cnt_sub <- count_curve(occ, sites = well_dated)
ps  <- phase_summary(fit, sd)
po  <- pairwise_overlap(occ, sd$site_names)

susq_save(plot_count_through_time(cnt, cnt_sub), file.path(figs, "fig1_count_through_time.png"), 7, 4.5)
susq_save(plot_occupation_intervals(ps, p$dates, p$curve), file.path(figs, "fig2_occupation_intervals.png"), 7, 8)
ord <- ps$site[order(ps$start_med)]; ord <- unique(ord[!is.na(ord)])
susq_save(plot_overlap_matrix(po, ord), file.path(figs, "fig3_overlap_matrix.png"), 7, 6.5)
susq_save(plot_duration_posteriors(fit, ps), file.path(figs, "fig4_duration_posteriors.png"), 7, 8)
susq_save(plot_calibrated_dates(p$dates, p$curve), file.path(figs, "fig5_calibrated_dates.png"), 7, 11)
cmp <- read.csv("output/susquehanna/model_comparison.csv", row.names = 1)
susq_save(plot_model_comparison(cmp), file.path(figs, "fig6_model_comparison.png"), 7, 3)
susq_save(plot_inbuilt_age(fit), file.path(figs, "fig7_inbuilt_age.png"), 7, 4)

write.csv(cnt, "output/susquehanna/count_curve.csv", row.names = FALSE)
write.csv(cnt_sub, "output/susquehanna/count_curve_well_dated.csv", row.names = FALSE)
write.csv(ps, "output/susquehanna/phase_summary.csv", row.names = FALSE)
write.csv(po, "output/susquehanna/pairwise_overlap.csv", row.names = FALSE)
cat("Figures written to", figs, "\n")
```

- [ ] **Step 5: Run it and look at every figure**

Run: `Rscript scripts/07_figures.R`
Then open each PNG. Check: no titles, Arial rendering, legible axis labels at 300 dpi, the reversal bands visible in figure 1, and figure 7 showing whether the inbuilt-age posteriors depart from the dashed uniform-prior line at 1/6.

- [ ] **Step 6: Commit**

```bash
git add R/susq_theme.R R/susq_plots.R scripts/07_figures.R tests/testthat/test_susq_theme.R output/susquehanna/figures output/susquehanna/*.csv
git commit -m "Add figure theme, seven result figures and the figure script"
```

---

### Task 11: Posterior predictive checks and the diagnostics figure

**Files:**
- Create: `scripts/06_posterior_predictive.R`
- Modify: `R/susq_plots.R` (append `plot_diagnostics`)

**Interfaces:**
- Produces `output/susquehanna/figures/fig8_diagnostics.png` and `output/susquehanna/ppc_summary.csv`.

The posterior predictive check compares each observed determination against replicates drawn from its own fitted mixture. Because the latent date was marginalised, replicates are generated in R by drawing a phase, then a variant, then a calendar date uniform in the phase, then a radiocarbon age from the calibration curve.

- [ ] **Step 1: Write the script**

```r
# scripts/06_posterior_predictive.R
for (f in c("R/susq_data.R","R/susq_calibration.R","R/susq_stan_data.R",
            "R/susq_fit.R","R/susq_metrics.R","R/susq_theme.R","R/susq_plots.R")) source(f)
set.seed(20260813)
p <- readRDS("output/susquehanna/prepared.rds")
fit <- readRDS("output/susquehanna/fits/multi.rds")
sd  <- readRDS("output/susquehanna/fits/multi_standata.rds")

st <- fit$draws("phase_start_calBP", format = "matrix")
en <- fit$draws("phase_end_calBP", format = "matrix")
z  <- fit$draws("z", format = "matrix")
D  <- nrow(st); n_rep <- 200L
draws <- sample(seq_len(D), n_rep)

yrep <- matrix(NA_real_, n_rep, sd$N)
for (r in seq_len(n_rep)) {
  d <- draws[r]
  for (i in seq_len(sd$N)) {
    ph <- z[d, i]
    t <- runif(1, en[d, ph], st[d, ph])
    mu <- approx(p$curve$calBP, p$curve$c14, xout = t, rule = 2)$y
    ce <- approx(p$curve$calBP, p$curve$err, xout = t, rule = 2)$y
    yrep[r, i] <- rnorm(1, mu, sqrt(p$dates$c14_error[i]^2 + ce^2))
  }
}

y <- p$dates$c14_age
ppc <- data.frame(
  det_id = p$dates$det_id, site = p$dates$site,
  material_class = p$dates$material_class,
  observed = y, rep_median = apply(yrep, 2, median),
  rep_lo = apply(yrep, 2, quantile, 0.025),
  rep_hi = apply(yrep, 2, quantile, 0.975))
ppc$covered <- ppc$observed >= ppc$rep_lo & ppc$observed <= ppc$rep_hi
write.csv(ppc, "output/susquehanna/ppc_summary.csv", row.names = FALSE)

cat(sprintf("PPC 95%% coverage overall: %.3f\n", mean(ppc$covered)))
print(tapply(ppc$covered, ppc$material_class, mean))
cat("\nCoverage far below 0.95 in one material class indicates that class's\n",
    "inbuilt-age or outlier treatment is not fitting.\n")

susq_save(plot_diagnostics(fit, ppc), "output/susquehanna/figures/fig8_diagnostics.png", 7, 6)
```

- [ ] **Step 2: Append the diagnostics plot to `R/susq_plots.R`**

```r
#' Figure 8. Supplementary diagnostics: posterior predictive coverage by
#' material class, and the distributions of R-hat and effective sample size.
plot_diagnostics <- function(fit, ppc) {
  s <- posterior::summarise_draws(fit$draws(c("mid", "dur", "mu_d", "sigma_d", "rho")),
                                  "rhat", "ess_bulk")
  a <- ggplot(ppc, aes(observed, rep_median, colour = material_class)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_linerange(aes(ymin = rep_lo, ymax = rep_hi), alpha = 0.35, linewidth = 0.3) +
    geom_point(size = 0.9) +
    scale_colour_manual(values = OKABE_ITO[c(4, 7, 3)], name = NULL) +
    labs(x = "Observed radiocarbon age (BP)",
         y = "Posterior predictive median (BP)") + theme_susq()
  b <- ggplot(s, aes(rhat)) + geom_histogram(bins = 40, fill = OKABE_ITO[6]) +
    geom_vline(xintercept = 1.01, linetype = "22", colour = OKABE_ITO[7]) +
    labs(x = "R-hat", y = "Parameters") + theme_susq()
  c_ <- ggplot(s, aes(ess_bulk)) + geom_histogram(bins = 40, fill = OKABE_ITO[6]) +
    geom_vline(xintercept = 400, linetype = "22", colour = OKABE_ITO[7]) +
    labs(x = "Bulk effective sample size", y = "Parameters") + theme_susq()
  patchwork::wrap_plots(a, patchwork::wrap_plots(b, c_, ncol = 2), ncol = 1,
                        heights = c(2, 1))
}
```

- [ ] **Step 3: Run it**

Run: `Rscript scripts/06_posterior_predictive.R`
Expected: overall coverage near 0.95. Coverage below about 0.85 in the `wood` class means the inbuilt-age model is not absorbing the offset, which must be investigated before the results are written up.

- [ ] **Step 4: Commit**

```bash
git add scripts/06_posterior_predictive.R R/susq_plots.R output/susquehanna/ppc_summary.csv output/susquehanna/figures/fig8_diagnostics.png
git commit -m "Add posterior predictive checks and diagnostics figure"
```

---

### Task 12: Results document and end-to-end run

**Files:**
- Create: `scripts/run_all.R`
- Create: `docs/susquehanna_results.md`

**Interfaces:**
- `scripts/run_all.R` runs the whole pipeline from a clean checkout.

- [ ] **Step 1: Write the runner**

```r
# scripts/run_all.R
#' Full pipeline. Expects setup_env.R to have been run once.
Sys.setenv(PATH = paste("/opt/homebrew/opt/gdal/bin", "/opt/homebrew/opt/geos/bin",
                        "/opt/homebrew/opt/proj/bin", Sys.getenv("PATH"), sep = ":"))
steps <- c("scripts/01_prepare_data.R",
           "scripts/02_prior_predictive.R",
           "scripts/03_simulation_recovery.R",
           "scripts/04_fit_models.R",
           "scripts/05_compare_models.R",
           "scripts/06_posterior_predictive.R",
           "scripts/07_figures.R")
for (s in steps) {
  cat("\n", strrep("=", 70), "\n", s, "\n", strrep("=", 70), "\n", sep = "")
  source(s, echo = FALSE)
}
cat("\nPipeline complete.\n")
```

- [ ] **Step 2: Run the tests end to end**

Run: `Rscript -e 'testthat::test_dir("tests/testthat")'`
Expected: all tests pass.

- [ ] **Step 3: Run the pipeline**

Run: `Rscript scripts/run_all.R 2>&1 | tee output/susquehanna/run_log.txt`

- [ ] **Step 4: Write `docs/susquehanna_results.md`**

Write the results in narrative prose, following the project writing conventions: no em dashes, no bullet lists in the prose, arguments leading and citations parenthetical, and no claim that outruns the evidence. The document must contain, in this order:

1. What was analysed: 164 determinations, 34 sites, how the merges and preprocessing were decided.
2. Whether the sampler and the simulation recovery passed, with the coverage numbers from `output/susquehanna/recovery_by_site_size.csv`. State plainly which site-size classes cannot support episode-count inference.
3. The contemporaneous-count result, read off `output/susquehanna/count_curve.csv`, reported both over all 34 sites and over the 22 well-dated ones, with the calibration reversals noted where they widen the interval.
4. Occupation durations, distinguishing sites where the posterior departs from the hierarchical prior from those where it does not.
5. Pairwise overlap, naming the specific pairs for which non-contemporaneity can be demonstrated, and stating how many pairs cannot be resolved.
6. Model comparison, with the caveat from spec section 5 that the sequential model was given its most favourable ordering.
7. Sensitivity to the duration prior across the 25, 60 and 150 year runs.
8. Limitations, including the `[CITE-CHECK: village occupation duration]` flag, the two-site basis for the wood inbuilt-age estimate, and the unresolved data-quality flags from `output/susquehanna/data_quality.csv`.

- [ ] **Step 5: Commit**

```bash
git add scripts/run_all.R docs/susquehanna_results.md output/susquehanna/run_log.txt
git commit -m "Add end-to-end runner and results document"
```

---

## Self-Review

**Spec coverage.** Section 2.2 site parsing is Task 1 and the existing `R/build_lookup_tables.R`. Section 2.3 material classes, Task 1. Section 2.4 preprocessing, Task 1. Section 2.5 quality flags, Task 1 and Task 12 step 4. Section 3.3 marginalisation, Tasks 2, 3 and 4. Section 3.4 inbuilt age, Tasks 2, 3 and 5, reported in figure 7. Section 3.5 outliers, Tasks 3 and 5. Section 3.6 full likelihood, Task 5. Section 3.7 duration prior and its sensitivity runs, Task 5 and Task 9 configs `multi_dur25` and `multi_dur150`. Section 4 derived quantities, Task 7. Section 5 model comparison, Tasks 6 and 9. Section 6 validation: prior predictive Task 9 step 2, simulation recovery Task 8, sampler diagnostics Task 5 `check_diagnostics`, posterior predictive Task 11, calibration cross-check against rcarbon Task 2, prior sensitivity Task 9. Section 7 figures, Tasks 10 and 11. Section 9 curve reversals, `CURVE_REVERSALS` in Task 10 and figure 1.

**Known gap.** Spec section 5 says exact 10-fold cross-validation replaces PSIS-LOO when more than five percent of Pareto-k values exceed 0.7. Task 9 step 4 detects and reports this but does not implement `scripts/05b_kfold.R`. Implement that script only if the check fires, since it costs ten refits per model and may not be needed.

**Type consistency.** `stan_data$site_names` is set in `build_stan_data` and consumed by `occupancy_from_arrays`, `pairwise_overlap` and `phase_summary`; it is stripped from the CmdStan payload in `fit_susq`. `phase_start` names a data array in Stan and a column prefix in the generated quantity `phase_start_calBP`; these are distinct and both are used consistently. `active`, `phase_start_calBP`, `phase_end_calBP` and `phase_duration` are produced by both `susq_phases.stan` and `susq_sequential.stan` with the same shapes, so `R/susq_metrics.R` works on either.
