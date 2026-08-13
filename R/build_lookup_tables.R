#' Generate the auditable lookup tables for site parsing and material class.
#'
#' Run once. The emitted CSVs in data/ are hand-audited afterwards and are the
#' authority used by the analysis. Re-running regenerates the proposal, so any
#' hand edits must be re-applied or the CSVs left alone.
#'
#' Decisions encoded here follow the design spec, section 2.2 and 2.3:
#'   docs/superpowers/specs/2026-08-13-susquehanna-contemporaneity-design.md

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
})

raw <- read_excel("data/Local Dates.xlsx")
names(raw)[1:6] <- c("site_raw", "material_raw", "lab_no",
                     "c14_age", "c14_error", "reference")
raw <- raw[, 1:8]
names(raw)[7:8] <- c("note_1", "note_2")

# ---------------------------------------------------------------------------
# Site parsing
# ---------------------------------------------------------------------------
# The site column concatenates site name and excavation context. Strip the
# context suffix, then apply the merges established from lab-number evidence
# and from the investigator's judgement (see spec section 2.2).

context_patterns <- c(
  "_F[ ._].*$", "_F\\d.*$", "_F$", "_PM.*$", "_post_mold.*$", "_Fea.*$",
  "_A horizon.*$", "_sample.*$", "_Sq.*$", "_sq.*$", "_from same.*$",
  "_shell.*$", "_fea.*$", "_burial.*$", "_Midden$", "_\\d+$", ",.*$",
  " \\(burial.*$", " Locus.*$"
)

strip_context <- function(x) {
  out <- str_trim(x)
  for (p in context_patterns) out <- str_remove(out, p)
  str_trim(out)
}

extract_context <- function(x) {
  base <- strip_context(x)
  ctx <- str_trim(str_remove(str_trim(x), fixed(base)))
  ctx <- str_remove(ctx, "^[_,]\\s*")
  ifelse(nchar(ctx) == 0, NA_character_, ctx)
}

# Merges. Left side is the stripped label, right side the analytical site.
#   Broome_Tech_PAF -> Broome_Tech: AA-31005/31006/31007 are consecutive lab
#     numbers spanning both labels, and Broome_Tech_1 cites Miroff_Kudrle_2017
#     along with every PAF date. PAF is the Public Archaeology Facility.
#   Chenango_Point_South -> Chenango_Point: Beta-35557/35558 are consecutive
#     and both trace to Wurst and Versaggi 1993. Merge directed by C. Lipo.
#   Sidney Hangar, Sidney AWOS -> Sidney Airport: three loci of one airport
#     complex investigated under separate projects. Merge directed by C. Lipo.
site_merges <- c(
  "Broome_Tech_PAF"      = "Broome_Tech",
  "Chenango_Point_South" = "Chenango_Point",
  "Sidney Hangar"        = "Sidney Airport",
  "Sidney AWOS"          = "Sidney Airport"
)

site_tbl <- raw %>%
  distinct(site_raw) %>%
  mutate(
    stripped = strip_context(site_raw),
    context  = extract_context(site_raw),
    site     = ifelse(stripped %in% names(site_merges),
                      site_merges[stripped], stripped),
    merged_from = ifelse(stripped %in% names(site_merges), stripped, NA_character_)
  ) %>%
  select(site_raw, site, context, merged_from) %>%
  arrange(site, context)

write.csv(site_tbl, "data/site_parse.csv", row.names = FALSE, na = "")

# ---------------------------------------------------------------------------
# Material classification
# ---------------------------------------------------------------------------
# Three classes. "short_lived" samples have no meaningful inbuilt age.
# "wood" samples may carry substantial inbuilt age. "indeterminate" covers
# descriptions that do not establish which, plus samples explicitly mixing
# wood and short-lived components; this class carries its own estimated
# inbuilt-age distribution rather than being assigned to wood by fiat.

short_lived <- c(
  "maize", "bean", "seed", "squash", "marshelder", "nutshell",
  "nutshell (juglans cinerea)", "butternut shell", "hazelnut shell",
  "bitternut shell", "charred buternut", "carbonized nutshell"
)

wood <- c(
  "probably wood", "wood", "wood charcoal", "charred wood", "wood (pm)",
  "oak wood charcoal",
  "wood charcoal, quercus sp.,ry1-5", "wood charcoal, quercus sp.,ry18-22",
  "wood charcoal, ulmus sp., ry1-5", "wood charcoal, ulmus sp., ry26-30"
)

indeterminate <- c(
  "charred material",          # description does not identify the taxon
  "carbonized material",       # as above
  "charcoal",                  # unqualified; probably wood but not stated
  "carbonized wood, nutshell, corn",  # mixed short-lived and wood
  "wood/nushell"                      # mixed short-lived and wood
)

classify <- function(m) {
  key <- str_squish(str_to_lower(m))
  case_when(
    key %in% short_lived   ~ "short_lived",
    key %in% wood          ~ "wood",
    key %in% indeterminate ~ "indeterminate",
    TRUE                   ~ NA_character_
  )
}

# Ring-year annotations: samples dated on counted tree rings, from the
# wiggle-matched sequences reported in Birch, Manning, Hart and Lorentzen 2024.
# ring_lo and ring_hi are the ring numbers spanned by the sample, counted
# outward from the innermost preserved ring.
ring_lo <- function(m) as.numeric(str_match(str_to_lower(m), "ry(\\d+)-\\d+")[, 2])
ring_hi <- function(m) as.numeric(str_match(str_to_lower(m), "ry\\d+-(\\d+)")[, 2])

mat_tbl <- raw %>%
  distinct(material_raw) %>%
  mutate(
    material_class = classify(material_raw),
    ring_lo        = ring_lo(material_raw),
    ring_hi        = ring_hi(material_raw),
    n_dates        = as.integer(table(str_squish(str_to_lower(raw$material_raw)))[
                       str_squish(str_to_lower(material_raw))])
  ) %>%
  arrange(material_class, desc(n_dates), material_raw)

stopifnot(!any(is.na(mat_tbl$material_class)))

write.csv(mat_tbl, "data/material_classes.csv", row.names = FALSE, na = "")

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
joined <- raw %>%
  left_join(site_tbl, by = "site_raw") %>%
  left_join(mat_tbl %>% select(material_raw, material_class, ring_lo, ring_hi),
            by = "material_raw")

cat("Determinations:", nrow(joined), "\n")
cat("Sites after merges:", n_distinct(joined$site), "\n\n")

cat("Material classes:\n")
print(table(joined$material_class))

cat("\nSites by number of determinations:\n")
n_by_site <- sort(table(joined$site), decreasing = TRUE)
print(n_by_site)

cat("\nPhase allocation (spec section 3.2):\n")
J <- ifelse(n_by_site >= 6, 4L, ifelse(n_by_site >= 3, 2L, 1L))
cat("  J=4 (n>=6):", sum(J == 4), "sites\n")
cat("  J=2 (n 3-5):", sum(J == 2), "sites\n")
cat("  J=1 (n<=2):", sum(J == 1), "sites\n")
cat("  total candidate phases:", sum(J), "\n")

cat("\nMaterial-class composition by site (identifies the inbuilt-age offsets):\n")
mix <- joined %>%
  count(site, material_class) %>%
  tidyr::pivot_wider(names_from = material_class, values_from = n, values_fill = 0) %>%
  filter(short_lived > 0, (wood > 0 | indeterminate > 0)) %>%
  arrange(desc(pmin(short_lived, wood + indeterminate)))
print(as.data.frame(mix))

cat("\n  sites informing the WOOD offset (short_lived > 0 & wood > 0):",
    sum(mix$short_lived > 0 & mix$wood > 0), "\n")
cat("  sites informing the INDETERMINATE offset (short_lived > 0 & indeterminate > 0):",
    sum(mix$short_lived > 0 & mix$indeterminate > 0), "\n")

cat("\nRing-year samples:\n")
print(as.data.frame(joined %>% filter(!is.na(ring_lo)) %>%
                    select(site, context, lab_no, c14_age, c14_error, ring_lo, ring_hi)))

cat("\nWrote data/site_parse.csv and data/material_classes.csv\n")
