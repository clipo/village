#' Load and preprocess the Susquehanna radiocarbon determinations.
#'
#' Applies the audited site and material lookup tables, pools the failed
#' Bates_66 replicate pair, and folds the two ring-counted wiggle-match pairs
#' into composite determinations. See spec sections 2.2 to 2.5.

suppressPackageStartupMessages({
  library(dplyr)
})

MATERIAL_CLASSES <- c("short_lived", "wood", "indeterminate")

#' Resolve a project-relative path regardless of caller's working directory.
#' testthat::test_file() runs tests with the working directory set to the
#' test file's own directory (tests/testthat), while analysis scripts run
#' from the project root; try both so the function's stated defaults work
#' in either context.
resolve_path <- function(path) {
  if (file.exists(path)) return(path)
  alt <- file.path("../..", path)
  if (file.exists(alt)) return(alt)
  stop("Cannot locate file: ", path, call. = FALSE)
}

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
  raw <- readxl::read_excel(resolve_path(xlsx))
  names(raw)[1:6] <- c("site_raw", "material_raw", "lab_no",
                       "c14_age", "c14_error", "reference")
  raw <- as.data.frame(raw[, 1:6])
  raw$c14_age <- as.numeric(raw$c14_age)
  raw$c14_error <- as.numeric(raw$c14_error)

  smap <- read.csv(resolve_path(site_map), stringsAsFactors = FALSE)
  mmap <- read.csv(resolve_path(material_map), stringsAsFactors = FALSE)

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
    # A duplicated base lab code would recycle the scalar pooled values across
    # unrelated rows, silently altering determinations this pipeline is
    # required to leave untouched. This dataset already contains one genuine
    # duplicate lab code (AA-41933), so fail loudly rather than assume.
    stopifnot(length(keep) == 1L)
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
