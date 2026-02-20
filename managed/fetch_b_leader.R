# fetch_b_leader.R -- Derive HLA-B -21 leader assignments from protein sequences
#
# Produces:
#   datasets/b_leader_assignments.csv
#
# Source: HLAtools protein alignment (position -21 of signal peptide)
#   M (methionine) = M-leader → efficient HLA-E loading → NKG2A-mediated NK education
#   T (threonine)  = T-leader → poor HLA-E loading → KIR-dependent NK education
#
# Also cross-validates against IPD Allele Query API (matching.b_leader field).
#
# Replaces family-level assignments from Petersdorf 2022 with allele-level
# sequence-derived classifications. Key differences from Petersdorf:
#   B*44: Petersdorf=M, sequence=T  |  B*45: Petersdorf=M, sequence=T
#   B*47: Petersdorf=M, sequence=T  |  B*56: Petersdorf=M, sequence=T
#   B*57: Petersdorf=M, sequence=T  |  B*58: Petersdorf=M, sequence=T

library(HLAtools)
library(httr)
library(jsonlite)
library(dplyr)
library(stringr)
library(readr)
library(here)

cat("=== fetch_b_leader.R ===\n")

# ============================================================
# 1. Derive -21 leader from HLAtools protein alignment
# ============================================================
cat("  Building HLA-B protein alignment from IMGT/HLA ...\n")

# Find the latest working IMGT/HLA version
B_align <- NULL
ipd_version <- NULL
for (v in c("3.63.0", "3.62.0", "3.61.0", "3.60.0", "3.58.0", "3.56.0")) {
  B_align <- tryCatch(buildAlignments("B", source = "AA", version = v),
                       error = function(e) NULL)
  if (!is.null(B_align)) {
    ipd_version <- B_align$B$Version
    cat(sprintf("  Alignment version: %s\n", ipd_version))
    break
  }
}

if (is.null(B_align)) stop("Could not build HLA-B protein alignment from HLAtools")

B_aa <- B_align$B$AA

# Extract position -21 and collapse to 2-field
b_leader <- B_aa %>%
  select(allele, trimmed_allele, pos_neg21 = `-21`) %>%
  filter(pos_neg21 %in% c("M", "T")) %>%
  mutate(allele_2field = trimmed_allele) %>%
  filter(!duplicated(allele_2field)) %>%
  transmute(
    allele = allele_2field,
    b_leader = pos_neg21,
    source = "HLAtools_protein_alignment",
    ipd_version = ipd_version,
    fetch_date = as.character(Sys.time())
  ) %>%
  arrange(allele)

cat(sprintf("  Derived %d unique 2-field alleles with -21 leader\n", nrow(b_leader)))
cat(sprintf("  Distribution: M=%d, T=%d\n",
            sum(b_leader$b_leader == "M"), sum(b_leader$b_leader == "T")))

# ============================================================
# 2. Cross-validate against IPD Allele Query API
# ============================================================
cat("  Cross-validating against IPD Allele Query API ...\n")

fetch_api_leaders <- function() {
  base_url <- "https://www.ebi.ac.uk/cgi-bin/ipd/api/allele"
  all_data <- list()
  page <- 1

  query_param <- 'startsWith(name,"B*")'
  url <- sprintf("%s?query=%s&limit=1000&fields=matching.b_leader",
                 base_url, URLencode(query_param, reserved = TRUE))

  repeat {
    if (page > 15) break  # safety limit
    resp <- tryCatch(GET(url), error = function(e) NULL)
    if (is.null(resp) || status_code(resp) != 200) break

    body <- content(resp, as = "text", encoding = "UTF-8")
    parsed <- fromJSON(body, flatten = TRUE)
    all_data <- c(all_data, list(parsed$data))

    next_path <- parsed$meta$`next`
    if (is.null(next_path) || is.na(next_path)) break
    url <- paste0(base_url, next_path)
    page <- page + 1
  }

  bind_rows(all_data)
}

api_raw <- tryCatch(fetch_api_leaders(), error = function(e) {
  cat(sprintf("  API unavailable: %s\n", e$message))
  data.frame()
})

if (nrow(api_raw) > 0) {
  api_leaders <- api_raw %>%
    filter(matching.b_leader %in% c("M", "T")) %>%
    mutate(allele_2field = str_extract(name, "B\\*\\d+:\\d+")) %>%
    filter(!is.na(allele_2field)) %>%
    group_by(allele_2field) %>%
    summarise(api_leader = names(sort(table(matching.b_leader), decreasing = TRUE))[1],
              .groups = "drop")

  cross_val <- b_leader %>%
    inner_join(api_leaders, by = c("allele" = "allele_2field")) %>%
    mutate(match = b_leader == api_leader)

  cat(sprintf("  Cross-validation: %d match, %d mismatch out of %d\n",
              sum(cross_val$match), sum(!cross_val$match), nrow(cross_val)))

  if (sum(!cross_val$match) > 0) {
    cat("  Mismatches (showing first 10):\n")
    cross_val %>% filter(!match) %>% select(allele, sequence = b_leader, api = api_leader) %>%
      head(10) %>% as.data.frame() %>% print(row.names = FALSE)
  }
} else {
  cat("  API cross-validation skipped (no data).\n")
}

# ============================================================
# 3. Write output
# ============================================================
output_path <- here("datasets", "b_leader_assignments.csv")
write_csv(b_leader, output_path)

cat(sprintf("\n  Wrote %s (%d rows)\n", basename(output_path), nrow(b_leader)))
cat(sprintf("  IPD-IMGT/HLA version: %s\n", ipd_version))
cat(sprintf("  Fetch date: %s\n", b_leader$fetch_date[1]))
