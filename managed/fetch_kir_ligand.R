# fetch_kir_ligand.R -- Fetch KIR ligand assignments from IPD-IMGT/HLA
#
# Produces:
#   datasets/bw4_bw6_classification.csv  (HLA-B: Bw4-80I, Bw4-80T, Bw6)
#   datasets/c1_c2_classification.csv     (HLA-C: C1, C2)
#
# Primary source: IPD-IMGT/HLA Allele Query API (bulk, paginated)
#   https://www.ebi.ac.uk/cgi-bin/ipd/api/allele?fields=matching.kir_ligand
#
# Cross-validation: HLAtools protein alignments (position 83 for Bw4/Bw6,
#   position 80 for C1/C2 and Bw4-80I/80T distinction)
#

library(httr)
library(jsonlite)
library(HLAtools)
library(dplyr)
library(stringr)
library(readr)
library(here)

cat("=== fetch_kir_ligand.R ===\n")

# ============================================================
# 1. Fetch from IPD Allele Query API (bulk)
# ============================================================
cat("  Fetching from IPD Allele Query API ...\n")

fetch_all_kir_ligand <- function(locus_prefix) {
  base_url <- "https://www.ebi.ac.uk/cgi-bin/ipd/api/allele"
  all_data <- list()
  page <- 1

  query_param <- sprintf('startsWith(name,"%s")', locus_prefix)
  url <- sprintf("%s?query=%s&limit=1000&fields=matching.kir_ligand",
                 base_url, URLencode(query_param, reserved = TRUE))

  repeat {
    cat(sprintf("    %s page %d ...\n", locus_prefix, page))
    resp <- GET(url)
    if (status_code(resp) != 200) {
      warning(sprintf("API returned status %d for %s page %d",
                       status_code(resp), locus_prefix, page))
      break
    }

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

collapse_to_2field <- function(raw_data, locus_letter) {
  pattern <- sprintf("%s\\*\\d+:\\d+", locus_letter)
  raw_data %>%
    filter(!is.na(matching.kir_ligand)) %>%
    mutate(allele_2field = str_extract(name, pattern)) %>%
    filter(!is.na(allele_2field)) %>%
    group_by(allele_2field) %>%
    summarise(
      kir_ligand = names(sort(table(matching.kir_ligand), decreasing = TRUE))[1],
      n_alleles_collapsed = n(),
      .groups = "drop"
    ) %>%
    mutate(source = "IPD_API") %>%
    arrange(allele_2field)
}

# Fetch HLA-B
hla_b_raw <- fetch_all_kir_ligand("B*")
hla_b_api <- collapse_to_2field(hla_b_raw, "B")
cat(sprintf("  HLA-B API: %d alleles fetched → %d unique 2-field\n",
            nrow(hla_b_raw), nrow(hla_b_api)))

# Fetch HLA-C
hla_c_raw <- fetch_all_kir_ligand("C*")
hla_c_api <- collapse_to_2field(hla_c_raw, "C")
cat(sprintf("  HLA-C API: %d alleles fetched → %d unique 2-field\n",
            nrow(hla_c_raw), nrow(hla_c_api)))

# ============================================================
# 2. Cross-validate with HLAtools protein alignments
# ============================================================
cat("  Cross-validating with HLAtools protein alignments ...\n")

# Find the latest working IMGT/HLA version
ipd_version <- NULL
for (v in c("3.63.0", "3.62.0", "3.61.0", "3.60.0", "3.58.0", "3.56.0")) {
  B_align <- tryCatch(buildAlignments("B", source = "AA", version = v),
                       error = function(e) NULL)
  if (!is.null(B_align)) {
    ipd_version <- B_align$B$Version
    cat(sprintf("  HLAtools alignment version: %s\n", ipd_version))
    break
  }
}

if (!is.null(B_align)) {
  B_aa <- B_align$B$AA
  # Derive Bw4/Bw6 from position 83 (R=Bw4, G=Bw6)
  b_seq <- B_aa %>%
    select(trimmed_allele, pos83 = `83`, pos80 = `80`) %>%
    filter(!duplicated(trimmed_allele), pos83 %in% c("R", "G")) %>%
    mutate(
      kir_ligand_seq = case_when(
        pos83 == "G" ~ "Bw6",
        pos83 == "R" & pos80 == "I" ~ "Bw4 - 80I",
        pos83 == "R" & pos80 == "T" ~ "Bw4 - 80T",
        pos83 == "R" ~ paste0("Bw4 - 80", pos80)  # rare variants
      ),
      allele_2field = trimmed_allele
    )

  # Cross-validate API vs sequence
  cross_val <- hla_b_api %>%
    inner_join(b_seq %>% select(allele_2field, kir_ligand_seq), by = "allele_2field") %>%
    mutate(match = kir_ligand == kir_ligand_seq)

  n_match <- sum(cross_val$match)
  n_mismatch <- sum(!cross_val$match)
  cat(sprintf("  Cross-validation (API vs sequence): %d match, %d mismatch out of %d\n",
              n_match, n_mismatch, nrow(cross_val)))

  if (n_mismatch > 0) {
    cat("  Mismatches (sequence is authoritative):\n")
    cross_val %>%
      filter(!match) %>%
      select(allele_2field, api = kir_ligand, sequence = kir_ligand_seq) %>%
      head(20) %>%
      as.data.frame() %>%
      print(row.names = FALSE)

    # Use sequence-derived classification where there's a mismatch
    hla_b_api <- hla_b_api %>%
      left_join(b_seq %>% select(allele_2field, kir_ligand_seq), by = "allele_2field") %>%
      mutate(
        kir_ligand = ifelse(!is.na(kir_ligand_seq), kir_ligand_seq, kir_ligand),
        source = ifelse(!is.na(kir_ligand_seq) & kir_ligand != kir_ligand_seq,
                        "sequence_corrected", source)
      ) %>%
      select(-kir_ligand_seq)
  }

  # Also do HLA-C cross-validation
  C_align <- tryCatch(buildAlignments("C", source = "AA", version = v),
                       error = function(e) NULL)
  if (!is.null(C_align)) {
    C_aa <- C_align$C$AA
    c_seq <- C_aa %>%
      select(trimmed_allele, pos80 = `80`) %>%
      filter(!duplicated(trimmed_allele), pos80 %in% c("N", "K")) %>%
      mutate(
        kir_ligand_seq = ifelse(pos80 == "N", "C1", "C2"),
        allele_2field = trimmed_allele
      )

    cross_val_c <- hla_c_api %>%
      inner_join(c_seq %>% select(allele_2field, kir_ligand_seq), by = "allele_2field") %>%
      mutate(match = kir_ligand == kir_ligand_seq)

    cat(sprintf("  HLA-C cross-validation: %d match, %d mismatch out of %d\n",
                sum(cross_val_c$match), sum(!cross_val_c$match), nrow(cross_val_c)))
  }
}

# ============================================================
# 3. Write output
# ============================================================
fetch_date <- as.character(Sys.time())

# ipd_version is resolved in section 2 from the HLAtools alignment and was
# previously only printed to stdout. It belongs in the tables themselves:
# downstream consumers pin a snapshot of these files and have to be able to
# state which IPD release produced the classifications they used.
#
# Assigned before fetch_date so the column order matches the other managed
# tables (... source, ipd_version, fetch_date). NULL would silently drop the
# column rather than write an empty one, so an unresolved version is recorded
# as NA -- absent, and visibly so.
ipd_version_out <- if (is.null(ipd_version)) NA_character_ else ipd_version
hla_b_api$ipd_version <- ipd_version_out
hla_c_api$ipd_version <- ipd_version_out

hla_b_api$fetch_date <- fetch_date
hla_c_api$fetch_date <- fetch_date

# `allele` is the uniform key across every managed table. `allele_2field` is
# kept as an exact duplicate for consumers written against the old name; it is
# deprecated and will be removed in v3.0.0.
#
# kir_ligand_code is the identifier form of kir_ligand -- no spaces or
# punctuation, so it is safe in a feature id, column name or model matrix.
# Computed here rather than in collapse_to_2field() because the cross-
# validation above can rewrite kir_ligand for sequence-corrected alleles, and
# the code must follow the corrected value.
add_key_and_code <- function(d) {
  d$allele <- d$allele_2field
  d$kir_ligand_code <- ifelse(
    grepl("^Bw4 - ", d$kir_ligand), sub("Bw4 - ", "Bw4_", d$kir_ligand), d$kir_ligand)
  front <- c("allele", "allele_2field", "kir_ligand", "kir_ligand_code")
  d[, c(front, setdiff(names(d), front))]
}

hla_b_api <- add_key_and_code(hla_b_api)
hla_c_api <- add_key_and_code(hla_c_api)

bw4_output <- here("datasets", "bw4_bw6_classification.csv")
c1c2_output <- here("datasets", "c1_c2_classification.csv")

write_csv(hla_b_api, bw4_output)
write_csv(hla_c_api, c1c2_output)

cat(sprintf("\n  Wrote %s (%d rows)\n", basename(bw4_output), nrow(hla_b_api)))
cat(sprintf("  Distribution: %s\n",
            paste(names(table(hla_b_api$kir_ligand)), table(hla_b_api$kir_ligand),
                  sep = "=", collapse = ", ")))
cat(sprintf("  Wrote %s (%d rows)\n", basename(c1c2_output), nrow(hla_c_api)))
cat(sprintf("  Distribution: %s\n",
            paste(names(table(hla_c_api$kir_ligand)), table(hla_c_api$kir_ligand),
                  sep = "=", collapse = ", ")))

cat(sprintf("\n  IPD-IMGT/HLA version: %s\n", ifelse(!is.null(ipd_version), ipd_version, "unknown")))
cat(sprintf("  Fetch date: %s\n", fetch_date))
