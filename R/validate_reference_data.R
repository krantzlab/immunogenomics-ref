# validate_reference_data.R -- Validation checks for reference_data module
#
# Run after any update to managed/ or curated/ data.
# Checks structural integrity, expected values, cross-references,
# and backward-compatible vector getter output.
#
# Usage (standalone):
#   source("R/validate_reference_data.R")
#   validate_all()   # returns TRUE if all checks pass
# Usage (as submodule):
#   source(here("reference_data", "R", "validate_reference_data.R"))
#   validate_all()

library(here)
library(readr)
library(dplyr)
library(stringr)

# Auto-detect: standalone repo or submodule
.validate_source_path <- if (file.exists(here::here("reference_data", "R", "load_reference_data.R"))) {
  here::here("reference_data", "R", "load_reference_data.R")
} else {
  here::here("R", "load_reference_data.R")
}
source(.validate_source_path)

# --- Individual validators ---

validate_bw4_bw6 <- function() {
  errors <- character()
  data <- tryCatch(load_bw4_classification(), error = function(e) {
    errors <<- c(errors, conditionMessage(e))
    tibble()
  })
  if (nrow(data) == 0) return(errors)

  # Known alleles must be present
  must_have <- c("B*44:02", "B*07:02", "B*27:05", "B*57:01")
  missing <- setdiff(must_have, data$allele_2field)
  if (length(missing) > 0) {
    errors <- c(errors, paste("Missing expected alleles:", paste(missing, collapse = ", ")))
  }

  # Spot-check known assignments
  b4402 <- data %>% filter(allele_2field == "B*44:02") %>% pull(kir_ligand)
  if (length(b4402) > 0 && b4402 != "Bw4 - 80T") {
    errors <- c(errors, sprintf("B*44:02 should be Bw4-80T, got: %s", b4402))
  }

  b5701 <- data %>% filter(allele_2field == "B*57:01") %>% pull(kir_ligand)
  if (length(b5701) > 0 && b5701 != "Bw4 - 80I") {
    errors <- c(errors, sprintf("B*57:01 should be Bw4-80I, got: %s", b5701))
  }

  errors
}

validate_bw4_80i <- function() {
  errors <- character()
  data <- tryCatch(load_bw4_80i_classification(), error = function(e) {
    errors <<- c(errors, conditionMessage(e))
    tibble()
  })
  if (nrow(data) == 0) return(errors)

  # Must have both loci
  if (!all(c("A", "B") %in% data$locus)) {
    errors <- c(errors, "Expected both HLA-A and HLA-B entries")
  }

  # Spot-check: B*57:01 must be Bw4-80I, B*44:02 must be Bw4-80T
  b5701 <- data %>% filter(allele_2field == "B*57:01")
  if (nrow(b5701) > 0 && b5701$kir_ligand != "Bw4 - 80I") {
    errors <- c(errors, sprintf("B*57:01 should be Bw4-80I, got: %s", b5701$kir_ligand))
  }

  b4402 <- data %>% filter(allele_2field == "B*44:02")
  if (nrow(b4402) > 0 && b4402$kir_ligand != "Bw4 - 80T") {
    errors <- c(errors, sprintf("B*44:02 should be Bw4-80T, got: %s", b4402$kir_ligand))
  }

  # A*32:01 must be Bw4-80I
  a3201 <- data %>% filter(allele_2field == "A*32:01")
  if (nrow(a3201) > 0 && a3201$kir_ligand != "Bw4 - 80I") {
    errors <- c(errors, sprintf("A*32:01 should be Bw4-80I, got: %s", a3201$kir_ligand))
  }

  errors
}

validate_c1_c2 <- function() {
  errors <- character()
  data <- tryCatch(load_c1_c2_classification(), error = function(e) {
    errors <<- c(errors, conditionMessage(e))
    tibble()
  })
  if (nrow(data) == 0) return(errors)

  # Spot-check
  c0501 <- data %>% filter(allele_2field == "C*05:01") %>% pull(kir_ligand)
  if (length(c0501) > 0 && c0501 != "C2") {
    errors <- c(errors, sprintf("C*05:01 should be C2, got: %s", c0501))
  }

  c0702 <- data %>% filter(allele_2field == "C*07:02") %>% pull(kir_ligand)
  if (length(c0702) > 0 && c0702 != "C1") {
    errors <- c(errors, sprintf("C*07:02 should be C1, got: %s", c0702))
  }

  errors
}

validate_b_leader <- function() {
  errors <- character()
  data <- tryCatch(load_b_leader_assignments(), error = function(e) {
    errors <<- c(errors, conditionMessage(e))
    tibble()
  })
  if (nrow(data) == 0) return(errors)

  # Spot-check against IPD-IMGT/HLA protein sequence data
  b4402 <- data %>% filter(allele == "B*44:02") %>% pull(b_leader)
  if (length(b4402) > 0 && b4402 != "T") {
    errors <- c(errors, sprintf("B*44:02 should be T-leader (IPD), got: %s", b4402))
  }

  b0702 <- data %>% filter(allele == "B*07:02") %>% pull(b_leader)
  if (length(b0702) > 0 && b0702 != "M") {
    errors <- c(errors, sprintf("B*07:02 should be M-leader, got: %s", b0702))
  }

  errors
}

validate_kir3dl1_expression <- function() {
  errors <- character()
  data <- tryCatch(load_kir3dl1_expression(), error = function(e) {
    errors <<- c(errors, conditionMessage(e))
    tibble()
  })
  if (nrow(data) == 0) return(errors)

  # *005 must be high (key classification decision)
  a005 <- data %>% filter(allele == "005") %>% pull(expression_class)
  if (length(a005) > 0 && a005 != "high") {
    errors <- c(errors, sprintf("*005 should be high expression, got: %s", a005))
  }

  # *004 must be low (not null)
  a004 <- data %>% filter(allele == "004") %>% pull(expression_class)
  if (length(a004) > 0 && a004 != "low") {
    errors <- c(errors, sprintf("*004 should be low expression, got: %s", a004))
  }

  # *054 must be null
  a054 <- data %>% filter(allele == "054") %>% pull(expression_class)
  if (length(a054) > 0 && a054 != "null") {
    errors <- c(errors, sprintf("*054 should be null expression, got: %s", a054))
  }

  errors
}

validate_kir3dl1_binding <- function() {
  errors <- character()
  data <- tryCatch(load_kir3dl1_binding(), error = function(e) {
    errors <<- c(errors, conditionMessage(e))
    tibble()
  })
  if (nrow(data) == 0) return(errors)

  # Must have all 9 KIR3DL1 allotypes
  expected_kir <- c("001", "002", "004", "005", "008", "009", "015", "020", "029")
  missing_kir <- setdiff(expected_kir, unique(data$kir3dl1_allele))
  if (length(missing_kir) > 0) {
    errors <- c(errors, sprintf("Missing KIR3DL1 allotypes: %s", paste(missing_kir, collapse = ", ")))
  }

  # Must have HLA alleles from all 3 loci
  has_a <- any(str_starts(data$hla_allele, "A\\*"))
  has_b <- any(str_starts(data$hla_allele, "B\\*"))
  has_c <- any(str_starts(data$hla_allele, "C\\*"))
  if (!has_a) errors <- c(errors, "No HLA-A alleles in binding data")
  if (!has_b) errors <- c(errors, "No HLA-B alleles in binding data")
  if (!has_c) errors <- c(errors, "No HLA-C alleles in binding data")

  # Spot-check: *008 vs B*57:01 should have highest MFI overall
  kir008_b5701 <- data %>% filter(kir3dl1_allele == "008", hla_allele == "B*57:01") %>% pull(mfi_mean)
  if (length(kir008_b5701) > 0) {
    max_mfi <- max(data$mfi_mean, na.rm = TRUE)
    if (kir008_b5701 != max_mfi) {
      errors <- c(errors, sprintf("*008 vs B*57:01 MFI (%.1f) should be highest (%.1f)", kir008_b5701, max_mfi))
    }
  }

  # Spot-check: A*32:01 must be present for multiple KIR allotypes
  a3201_data <- data %>% filter(hla_allele == "A*32:01")
  if (nrow(a3201_data) < 5) {
    errors <- c(errors, sprintf("Expected >= 5 KIR allotypes measured against A*32:01, got %d", nrow(a3201_data)))
  }

  # Total should be >= 300 measurements
  if (nrow(data) < 300) {
    errors <- c(errors, sprintf("Expected >= 300 binding measurements, got %d", nrow(data)))
  }

  errors
}

validate_erap1_activity <- function() {
  errors <- character()
  data <- tryCatch(load_erap1_allotypes(), error = function(e) {
    errors <<- c(errors, conditionMessage(e))
    tibble()
  })
  if (nrow(data) == 0) return(errors)

  # Allotype 10 must be very_low activity
  a10 <- data %>% filter(allotype == 10) %>% pull(activity_class)
  if (length(a10) > 0 && a10 != "very_low") {
    errors <- c(errors, sprintf("Allotype 10 should be very_low, got: %s", a10))
  }

  # Allotype 2 must be high activity (most common globally)
  a2 <- data %>% filter(allotype == 2) %>% pull(activity_class)
  if (length(a2) > 0 && a2 != "high") {
    errors <- c(errors, sprintf("Allotype 2 should be high, got: %s", a2))
  }

  # Must have SNP haplotype columns
  snp_cols <- c("rs30187_pos528_ancA_derG", "rs27044_pos730_ancC_derG")
  missing_snp <- setdiff(snp_cols, names(data))
  if (length(missing_snp) > 0) {
    errors <- c(errors, paste("Missing SNP columns:", paste(missing_snp, collapse = ", ")))
  }

  errors
}

validate_erap2_expression <- function() {
  errors <- character()
  data <- tryCatch(load_erap2_haplotypes(), error = function(e) {
    errors <<- c(errors, conditionMessage(e))
    tibble()
  })
  if (nrow(data) == 0) return(errors)

  if (nrow(data) != 3) {
    errors <- c(errors, sprintf("Expected exactly 3 haplotypes (HapA/HapB/HapC), got %d", nrow(data)))
  }

  # HapB (G allele) must be non-expressing
  hapb <- data %>% filter(haplotype == "HapB") %>% pull(expression_status)
  if (length(hapb) > 0 && hapb != "non_expressing") {
    errors <- c(errors, sprintf("HapB should be non_expressing, got: %s", hapb))
  }

  # HapA (A allele) must be expressing
  hapa <- data %>% filter(haplotype == "HapA") %>% pull(expression_status)
  if (length(hapa) > 0 && hapa != "expressing") {
    errors <- c(errors, sprintf("HapA should be expressing, got: %s", hapa))
  }

  errors
}

validate_hla_a_expression <- function() {
  errors <- character()
  data <- tryCatch(load_hla_a_expression(), error = function(e) {
    errors <<- c(errors, conditionMessage(e))
    tibble()
  })
  if (nrow(data) == 0) return(errors)

  # Must have >= 20 lineages
  if (nrow(data) < 20) {
    errors <- c(errors, sprintf("Expected >= 20 HLA-A lineages, got %d", nrow(data)))
  }

  # Common lineages must be present
  must_have <- c("A*01", "A*02", "A*03", "A*24")
  missing <- setdiff(must_have, data$hla_a_lineage)
  if (length(missing) > 0) {
    errors <- c(errors, paste("Missing expected lineages:", paste(missing, collapse = ", ")))
  }

  # Z-scores should be roughly centered (mean near 0)
  mean_z <- mean(data$expression_z_score)
  if (abs(mean_z) > 0.5) {
    errors <- c(errors, sprintf("Mean z-score = %.3f (expected near 0)", mean_z))
  }

  # Spot-check: A*24 should be high, A*74 should be low
  a24 <- data %>% filter(hla_a_lineage == "A*24") %>% pull(expression_class)
  if (length(a24) > 0 && a24 != "high") {
    errors <- c(errors, sprintf("A*24 should be high expression, got: %s", a24))
  }

  a74 <- data %>% filter(hla_a_lineage == "A*74") %>% pull(expression_class)
  if (length(a74) > 0 && a74 != "low") {
    errors <- c(errors, sprintf("A*74 should be low expression, got: %s", a74))
  }

  errors
}

validate_hla_c_expression <- function() {
  errors <- character()
  data <- tryCatch(load_hla_c_expression(), error = function(e) {
    errors <<- c(errors, conditionMessage(e))
    tibble()
  })
  if (nrow(data) == 0) return(errors)

  # Must have >= 14 allotypes
  if (nrow(data) < 14) {
    errors <- c(errors, sprintf("Expected >= 14 HLA-C allotypes, got %d", nrow(data)))
  }

  # Common allotypes must be present
  must_have <- c("C*04", "C*06", "C*07", "C*14")
  missing <- setdiff(must_have, data$hla_c_allotype)
  if (length(missing) > 0) {
    errors <- c(errors, paste("Missing expected allotypes:", paste(missing, collapse = ", ")))
  }

  # MFI values must be positive
  if (any(data$expression_mfi <= 0, na.rm = TRUE)) {
    errors <- c(errors, "MFI values must be positive")
  }

  # Spot-check: C*07 should be low (lowest MFI), C*14 should be high (highest MFI)
  c07 <- data %>% filter(hla_c_allotype == "C*07") %>% pull(expression_class)
  if (length(c07) > 0 && c07 != "low") {
    errors <- c(errors, sprintf("C*07 should be low expression, got: %s", c07))
  }

  c14 <- data %>% filter(hla_c_allotype == "C*14") %>% pull(expression_class)
  if (length(c14) > 0 && c14 != "high") {
    errors <- c(errors, sprintf("C*14 should be high expression, got: %s", c14))
  }

  # C*14 should have the highest MFI
  c14_mfi <- data %>% filter(hla_c_allotype == "C*14") %>% pull(expression_mfi)
  max_mfi <- max(data$expression_mfi, na.rm = TRUE)
  if (length(c14_mfi) > 0 && c14_mfi != max_mfi) {
    errors <- c(errors, sprintf("C*14 MFI (%s) should be highest (%s)", c14_mfi, max_mfi))
  }

  errors
}

validate_hla_divergence <- function() {
  errors <- character()
  data <- tryCatch(load_hla_divergence(), error = function(e) {
    errors <<- c(errors, conditionMessage(e))
    tibble()
  })
  if (nrow(data) == 0) return(errors)

  # Must have all 3 loci
  missing_loci <- setdiff(c("A", "B", "C"), unique(data$locus))
  if (length(missing_loci) > 0) {
    errors <- c(errors, sprintf("Missing loci: %s", paste(missing_loci, collapse = ", ")))
  }

  # FD scores must be in [0, 1]
  if (any(data$fd_score < 0 | data$fd_score > 1, na.rm = TRUE)) {
    errors <- c(errors, "fd_score values outside [0, 1] range")
  }

  # Self-pairs must have fd_score == 0
  self_pairs <- data %>% filter(allele_1 == allele_2)
  if (nrow(self_pairs) > 0 && any(self_pairs$fd_score != 0)) {
    errors <- c(errors, "Self-pairs (allele_1 == allele_2) should have fd_score = 0")
  }

  # Spot-check: A*01:01 vs A*01:01 = 0
  a0101 <- data %>% filter(allele_1 == "A*01:01", allele_2 == "A*01:01")
  if (nrow(a0101) > 0 && a0101$fd_score[1] != 0) {
    errors <- c(errors, sprintf("A*01:01 self-pair should be 0, got: %s", a0101$fd_score[1]))
  }

  errors
}

validate_erap_crosswalk <- function() {
  errors <- character()
  cw <- tryCatch(load_erap_crosswalk(), error = function(e) {
    errors <<- c(errors, conditionMessage(e))
    tibble()
  })
  if (nrow(cw) == 0) return(errors)

  activity <- tryCatch(load_erap1_allotypes(), error = function(e) {
    errors <<- c(errors, paste("load_erap1_allotypes():", conditionMessage(e)))
    tibble()
  })
  if (nrow(activity) == 0) return(errors)

  # Parse the activity-table SNP columns (e.g. rs2287987_pos349_ancA_derG) into
  # rsid -> (aa_position, ancestral residue, derived residue). The ancestral
  # residue is the value carried by allotype 1 (the ancestral allotype); the
  # derived residue is the other allele observed in that column.
  snp_cols <- grep("^rs", names(activity), value = TRUE)
  enc <- lapply(snp_cols, function(col) {
    rsid <- sub("_.*$", "", col)
    pos  <- as.integer(str_match(col, "_pos(\\d+)_")[, 2])
    vals <- unique(activity[[col]])
    anc  <- activity[[col]][activity$allotype == 1][1]
    der  <- setdiff(vals, anc)
    list(rsid = rsid, col = col, pos = pos, anc = anc, der = der)
  })
  names(enc) <- vapply(enc, function(e) e$rsid, character(1))

  # (a) rsID set: every ERAP1 activity SNP + rs2248374 must be present; the only
  #     extra row allowed is the optional rs2549782 (ERAP2 N392K candidate).
  required <- c(names(enc), "rs2248374")
  optional_allowed <- "rs2549782"
  missing_rs <- setdiff(required, cw$rsid)
  if (length(missing_rs) > 0) {
    errors <- c(errors, sprintf("crosswalk missing required rsIDs: %s", paste(missing_rs, collapse = ", ")))
  }
  extra_rs <- setdiff(cw$rsid, c(required, optional_allowed))
  if (length(extra_rs) > 0) {
    errors <- c(errors, sprintf("crosswalk has unexpected rsIDs: %s", paste(extra_rs, collapse = ", ")))
  }

  # (b) Each ERAP1 SNP's aa_position / ancestral / derived must match the
  #     activity-table encoding exactly (position and direction).
  for (e in enc) {
    row <- cw %>% filter(rsid == e$rsid)
    if (nrow(row) == 0) next  # already reported by (a)
    if (length(e$der) != 1) {
      errors <- c(errors, sprintf("%s: activity column %s does not have exactly 2 residues", e$rsid, e$col))
      next
    }
    if (is.na(row$aa_position) || row$aa_position != e$pos) {
      errors <- c(errors, sprintf("%s: aa_position %s != activity-table %s", e$rsid, row$aa_position, e$pos))
    }
    if (row$aa_ancestral != e$anc || row$aa_derived != e$der) {
      errors <- c(errors, sprintf("%s: crosswalk anc/der %s/%s != activity-table %s/%s",
                                  e$rsid, row$aa_ancestral, row$aa_derived, e$anc, e$der))
    }
  }

  # (c) ref_aa/alt_aa must be strand-consistent: the {ref_aa, alt_aa} pair must
  #     equal the {ancestral, derived} pair, and coding_strand must match the
  #     gene's orientation at 5q15 (ERAP1 minus, ERAP2 plus).
  gene_strand <- c(ERAP1 = "minus", ERAP2 = "plus")
  for (i in seq_len(nrow(cw))) {
    row <- cw[i, ]
    expect_strand <- gene_strand[[row$gene]]
    if (!is.null(expect_strand) && row$coding_strand != expect_strand) {
      errors <- c(errors, sprintf("%s (%s): coding_strand %s != expected %s",
                                  row$rsid, row$gene, row$coding_strand, expect_strand))
    }
    # Coding rows only (splice/non-coding rows carry NA amino acids).
    if (!is.na(row$ref_aa) && !is.na(row$alt_aa) &&
        !is.na(row$aa_ancestral) && !is.na(row$aa_derived)) {
      if (!setequal(c(row$ref_aa, row$alt_aa), c(row$aa_ancestral, row$aa_derived))) {
        errors <- c(errors, sprintf("%s: {ref_aa,alt_aa}={%s,%s} inconsistent with {anc,der}={%s,%s}",
                                    row$rsid, row$ref_aa, row$alt_aa, row$aa_ancestral, row$aa_derived))
      }
    }
  }

  # (d) No missing GRCh38 coordinates on any row.
  if (any(is.na(cw$grch38_chrom)) || any(is.na(cw$grch38_pos))) {
    bad <- cw$rsid[is.na(cw$grch38_chrom) | is.na(cw$grch38_pos)]
    errors <- c(errors, sprintf("crosswalk rows missing GRCh38 coordinates: %s", paste(bad, collapse = ", ")))
  }

  errors
}

# --- Backward-compatible getter validation ---

validate_getters <- function() {
  errors <- character()

  # Bw4 vectors
  bw4_80i <- tryCatch(get_bw4_80i_alleles(), error = function(e) {
    errors <<- c(errors, paste("get_bw4_80i_alleles():", conditionMessage(e)))
    character()
  })
  if (length(bw4_80i) > 0) {
    if (!"B*57:01" %in% bw4_80i) errors <- c(errors, "get_bw4_80i_alleles: missing B*57:01")
    if ("B*44:02" %in% bw4_80i)  errors <- c(errors, "get_bw4_80i_alleles: B*44:02 should NOT be 80I")
    if (!all(str_starts(bw4_80i, "B\\*"))) errors <- c(errors, "get_bw4_80i_alleles: non-B* alleles found")
  }

  bw4_80t <- tryCatch(get_bw4_80t_alleles(), error = function(e) {
    errors <<- c(errors, paste("get_bw4_80t_alleles():", conditionMessage(e)))
    character()
  })
  if (length(bw4_80t) > 0) {
    if (!"B*44:02" %in% bw4_80t) errors <- c(errors, "get_bw4_80t_alleles: missing B*44:02")
    if (!"B*27:05" %in% bw4_80t) errors <- c(errors, "get_bw4_80t_alleles: missing B*27:05 (should be 80T)")
  }

  a_bw4 <- tryCatch(get_hla_a_bw4_alleles(), error = function(e) {
    errors <<- c(errors, paste("get_hla_a_bw4_alleles():", conditionMessage(e)))
    character()
  })
  if (length(a_bw4) > 0) {
    if (!"A*32:01" %in% a_bw4) errors <- c(errors, "get_hla_a_bw4_alleles: missing A*32:01")
    if (!all(str_starts(a_bw4, "A\\*"))) errors <- c(errors, "get_hla_a_bw4_alleles: non-A* alleles found")
  }

  # No overlap between 80I and 80T
  overlap <- intersect(bw4_80i, bw4_80t)
  if (length(overlap) > 0) {
    errors <- c(errors, sprintf("Bw4-80I and 80T vectors overlap: %s",
                                paste(head(overlap, 5), collapse = ", ")))
  }

  # C1/C2 vectors
  c1 <- tryCatch(get_c1_alleles(), error = function(e) {
    errors <<- c(errors, paste("get_c1_alleles():", conditionMessage(e)))
    character()
  })
  c2 <- tryCatch(get_c2_alleles(), error = function(e) {
    errors <<- c(errors, paste("get_c2_alleles():", conditionMessage(e)))
    character()
  })
  if (length(c1) > 0 && length(c2) > 0) {
    c_overlap <- intersect(c1, c2)
    if (length(c_overlap) > 0) {
      errors <- c(errors, sprintf("C1 and C2 vectors overlap: %s",
                                  paste(head(c_overlap, 5), collapse = ", ")))
    }
  }

  # KIR3DL1 expression vectors: must be character, no overlap
  kir_high <- tryCatch(get_kir3dl1_high(), error = function(e) {
    errors <<- c(errors, paste("get_kir3dl1_high():", conditionMessage(e)))
    character()
  })
  kir_low <- tryCatch(get_kir3dl1_low(), error = function(e) {
    errors <<- c(errors, paste("get_kir3dl1_low():", conditionMessage(e)))
    character()
  })
  kir_null <- tryCatch(get_kir3dl1_null(), error = function(e) {
    errors <<- c(errors, paste("get_kir3dl1_null():", conditionMessage(e)))
    character()
  })
  if (length(kir_high) > 0) {
    if (!"005" %in% kir_high) errors <- c(errors, "get_kir3dl1_high: missing 005")
    if (!"001" %in% kir_high) errors <- c(errors, "get_kir3dl1_high: missing 001")
  }
  kir_all <- c(kir_high, kir_low, kir_null)
  if (length(kir_all) != length(unique(kir_all))) {
    errors <- c(errors, "KIR3DL1 expression vectors have duplicate alleles across tiers")
  }

  # KIR3DL1 MFI: must be named numeric
  mfi <- tryCatch(get_kir3dl1_mfi(), error = function(e) {
    errors <<- c(errors, paste("get_kir3dl1_mfi():", conditionMessage(e)))
    numeric()
  })
  if (length(mfi) > 0) {
    if (!is.numeric(mfi)) errors <- c(errors, "get_kir3dl1_mfi: not numeric")
    if (is.null(names(mfi))) errors <- c(errors, "get_kir3dl1_mfi: not named")
    if (!"008" %in% names(mfi)) errors <- c(errors, "get_kir3dl1_mfi: missing 008")
  }

  # Cross-validation: every KIR3DL1 allotype in binding must also be in expression
  binding_kir <- tryCatch(
    unique(load_kir3dl1_binding()$kir3dl1_allele), error = function(e) character()
  )
  expr_alleles <- tryCatch(
    load_kir3dl1_expression()$allele, error = function(e) character()
  )
  if (length(binding_kir) > 0 && length(expr_alleles) > 0) {
    missing_expr <- setdiff(binding_kir, expr_alleles)
    if (length(missing_expr) > 0) {
      errors <- c(errors, sprintf("KIR3DL1 allotypes in binding but not expression: %s",
                                  paste(missing_expr, collapse = ", ")))
    }
  }

  # B leader lookup: must be named character
  bl <- tryCatch(get_b_leader_lookup(), error = function(e) {
    errors <<- c(errors, paste("get_b_leader_lookup():", conditionMessage(e)))
    character()
  })
  if (length(bl) > 0) {
    if (!is.character(bl)) errors <- c(errors, "get_b_leader_lookup: not character")
    if (is.null(names(bl))) errors <- c(errors, "get_b_leader_lookup: not named")
    if ("B*07:02" %in% names(bl) && bl["B*07:02"] != "M") {
      errors <- c(errors, sprintf("get_b_leader_lookup: B*07:02 should be M, got %s", bl["B*07:02"]))
    }
  }

  errors
}

# --- Master validator ---

#' Run all validation checks
#' @param verbose logical; if TRUE, print results to console
#' @return TRUE if all pass, FALSE otherwise
validate_all <- function(verbose = TRUE) {
  checks <- list(
    bw4_bw6             = validate_bw4_bw6,
    bw4_80i             = validate_bw4_80i,
    c1_c2               = validate_c1_c2,
    b_leader            = validate_b_leader,
    kir3dl1_expression  = validate_kir3dl1_expression,
    kir3dl1_binding     = validate_kir3dl1_binding,
    erap1_activity      = validate_erap1_activity,
    erap2_expression    = validate_erap2_expression,
    erap_crosswalk      = validate_erap_crosswalk,
    hla_a_expression    = validate_hla_a_expression,
    hla_c_expression    = validate_hla_c_expression,
    hla_divergence      = validate_hla_divergence,
    backward_compat     = validate_getters
  )

  all_pass <- TRUE
  for (name in names(checks)) {
    errors <- checks[[name]]()
    if (length(errors) == 0) {
      if (verbose) cat(sprintf("  PASS: %s\n", name))
    } else {
      all_pass <- FALSE
      if (verbose) {
        cat(sprintf("  FAIL: %s\n", name))
        for (e in errors) cat(sprintf("    - %s\n", e))
      }
    }
  }

  if (verbose) {
    if (all_pass) {
      cat("\n  All reference data checks passed.\n")
    } else {
      cat("\n  Some checks FAILED. Review errors above.\n")
    }
  }

  invisible(all_pass)
}
