# load_reference_data.R -- Loader functions for the reference_data module
#
# Two kinds of functions:
#   1. load_*() -- Read from datasets/ with column and value validation.
#                  Return tidy tibbles for flexible downstream use.
#   2. get_*()  -- Backward-compatible vector extractors that return the same
#                  formats as hardcoded constants in downstream projects.
#
# Usage (standalone):
#   source("R/load_reference_data.R")
# Usage (as submodule):
#   source(here("reference_data", "R", "load_reference_data.R"))
#
#   # Replace hardcoded constants:
#   BW4_80I_ALLELES        <- get_bw4_80i_alleles()
#   BW4_80T_ALLELES        <- get_bw4_80t_alleles()
#   HLA_A_BW4_80I_ALLELES  <- get_hla_a_bw4_alleles()
#   C1_ALLELES             <- get_c1_alleles()
#   C2_ALLELES             <- get_c2_alleles()
#   KIR3DL1_HIGH           <- get_kir3dl1_high()
#   KIR3DL1_LOW            <- get_kir3dl1_low()
#   KIR3DL1_NULL           <- get_kir3dl1_null()
#   KIR3DL1_A3201_MFI      <- get_kir3dl1_mfi()
#   KIR3DL1_BINDING_STRONG <- get_kir3dl1_binding_strong()
#   ... etc.

library(readr)
library(dplyr)
library(stringr)
library(here)

# ============================================================
# Internal helpers
# ============================================================

.refdata_root <- function() {
  # Context 2: submodule inside parent project
  submod <- here::here("reference_data", "datasets")
  if (dir.exists(submod)) return(here::here("reference_data"))
  # Context 1: standalone repo
  here::here()
}

.refdata_path <- function(filename) {
  file.path(.refdata_root(), "datasets", filename)
}

.load_csv <- function(filename, required_cols, context) {
  path <- .refdata_path(filename)
  if (!file.exists(path)) {
    stop(sprintf("%s not found at %s\n  Run the appropriate managed/ script first.", filename, path),
         call. = FALSE)
  }
  data <- read_csv(path, show_col_types = FALSE)
  missing <- setdiff(required_cols, names(data))
  if (length(missing) > 0) {
    stop(sprintf("[%s] Missing required columns: %s", context, paste(missing, collapse = ", ")),
         call. = FALSE)
  }
  data
}

# ============================================================
# Track 1: Managed (from authoritative databases)
# ============================================================

#' Load Bw4/Bw6 classification for HLA-B alleles
#' @return tibble: allele_2field, kir_ligand, n_alleles_collapsed, source, ipd_version, fetch_date
load_bw4_classification <- function() {
  data <- .load_csv("bw4_bw6_classification.csv",
                     c("allele_2field", "kir_ligand"),
                     "bw4_classification")
  valid <- c("Bw4 - 80I", "Bw4 - 80T", "Bw4 - 80N", "Bw6")
  bad <- setdiff(unique(data$kir_ligand), valid)
  if (length(bad) > 0) {
    stop(sprintf("[bw4_classification] Unexpected kir_ligand values: %s", paste(bad, collapse = ", ")),
         call. = FALSE)
  }
  stopifnot("All alleles must start with B*" = all(str_starts(data$allele_2field, "B\\*")))
  data
}

#' Load Bw4 position-80 classification (HLA-B + HLA-A)
#'
#' Null alleles are excluded by design: they are not expressed, so they present
#' no epitope and cannot be KIR ligands. As a consequence every row carries a
#' readable position 83, and `classification_status` has exactly one non-
#' classified value -- `non_canonical_83`, meaning the allele carries neither
#' canonical epitope. There is no "sequence unavailable" state.
#'
#' Coverage differs by locus. HLA-B is complete: every expressed 2-field allele
#' is present. HLA-A is partial by construction -- only Bw4 carriers are
#' listed -- so an HLA-A allele's absence means "not Bw4" and is informative,
#' while an HLA-B allele's absence means the snapshot predates it.
#'
#' @return tibble: allele_2field, locus, kir_ligand, kir_ligand_code,
#'   classification_status, pos80, bw4_motif_77_83, api_kir_ligand, source,
#'   ipd_version, fetch_date
load_bw4_80i_classification <- function() {
  data <- .load_csv("bw4_80i_classification.csv",
                     c("allele_2field", "locus", "kir_ligand", "kir_ligand_code",
                       "classification_status", "pos80"),
                     "bw4_80i_classification")
  stopifnot("locus must be A or B" = all(data$locus %in% c("A", "B")))

  bad_status <- setdiff(unique(data$classification_status),
                        c("classified", "non_canonical_83"))
  if (length(bad_status) > 0) {
    stop(sprintf("[bw4_80i_classification] Unexpected classification_status: %s",
                 paste(bad_status, collapse = ", ")), call. = FALSE)
  }

  # The identifier column is what downstream code builds feature names from, so
  # a space or punctuation slipping in is a downstream break, not a cosmetic one.
  bad_code <- unique(data$kir_ligand_code[!str_detect(data$kir_ligand_code, "^[A-Za-z0-9_]+$")])
  if (length(bad_code) > 0) {
    stop(sprintf("[bw4_80i_classification] kir_ligand_code must be syntax-safe, got: %s",
                 paste(bad_code, collapse = ", ")), call. = FALSE)
  }
  data
}

#' Load C1/C2 classification for HLA-C alleles
#' @return tibble: allele_2field, kir_ligand, n_alleles_collapsed, source, ipd_version, fetch_date
load_c1_c2_classification <- function() {
  data <- .load_csv("c1_c2_classification.csv",
                     c("allele_2field", "kir_ligand"),
                     "c1_c2_classification")
  bad <- setdiff(unique(data$kir_ligand), c("C1", "C2"))
  if (length(bad) > 0) {
    stop(sprintf("[c1_c2_classification] Unexpected kir_ligand values: %s", paste(bad, collapse = ", ")),
         call. = FALSE)
  }
  data
}

#' Load HLA-B -21 leader peptide (M/T) assignments
#' @return tibble: allele, b_leader, source, ipd_version, fetch_date
load_b_leader_assignments <- function() {
  data <- .load_csv("b_leader_assignments.csv",
                     c("allele", "b_leader"),
                     "b_leader_assignments")
  bad <- setdiff(unique(data$b_leader), c("M", "T"))
  if (length(bad) > 0) {
    stop(sprintf("[b_leader_assignments] Unexpected b_leader values: %s", paste(bad, collapse = ", ")),
         call. = FALSE)
  }
  data
}

# ============================================================
# Track 2: Curated (from literature)
# ============================================================

#' Load KIR3DL1 expression classification
#' @return tibble: allele, expression_class, confidence_tier, source_doi, source_detail, evidence_note
load_kir3dl1_expression <- function() {
  data <- .load_csv("kir3dl1_expression.csv",
                     c("allele", "expression_class", "confidence_tier", "source_doi"),
                     "kir3dl1_expression")
  bad <- setdiff(unique(data$expression_class), c("high", "low", "null"))
  if (length(bad) > 0) {
    stop(sprintf("[kir3dl1_expression] Unexpected expression_class: %s", paste(bad, collapse = ", ")),
         call. = FALSE)
  }
  data
}

#' Load KIR3DL1-HLA binding data (full matrix from bead-based assay)
#' @param kir3dl1 optional; filter to specific KIR3DL1 allotype (e.g., "001")
#' @param hla optional; filter to specific HLA allele (e.g., "A*32:01")
#' @return tibble: kir3dl1_allele, hla_allele, mfi_mean, confidence_tier, source_doi, source_detail
load_kir3dl1_binding <- function(kir3dl1 = NULL, hla = NULL) {
  data <- .load_csv("kir3dl1_hla_binding.csv",
                     c("kir3dl1_allele", "hla_allele", "mfi_mean", "confidence_tier", "source_doi"),
                     "kir3dl1_binding")
  stopifnot("mfi_mean must be numeric" = is.numeric(data$mfi_mean))
  if (!is.null(kir3dl1)) {
    data <- data %>% filter(kir3dl1_allele %in% kir3dl1)
  }
  if (!is.null(hla)) {
    data <- data %>% filter(hla_allele %in% hla)
  }
  data
}

#' Load ERAP1 allotype activity with SNP haplotypes
#' @return tibble: allotype, activity_class, rs*_pos*_anc*_der* columns, freq_EUR, confidence_tier, source_doi, evidence_note
load_erap1_allotypes <- function() {
  data <- .load_csv("erap1_allotype_activity.csv",
                     c("allotype", "activity_class", "freq_EUR", "confidence_tier", "source_doi"),
                     "erap1_allotypes")
  valid_classes <- c("high", "moderate", "low_to_moderate", "low", "very_low")
  bad <- setdiff(unique(data$activity_class), valid_classes)
  if (length(bad) > 0) {
    stop(sprintf("[erap1_allotypes] Unexpected activity_class: %s", paste(bad, collapse = ", ")),
         call. = FALSE)
  }
  stopifnot("Expected >= 10 ERAP1 allotypes" = nrow(data) >= 10)
  data
}

#' Load ERAP2 haplotype expression
#' @return tibble: haplotype, tag_snp_rs2248374, expression_status, expression_level,
#'   functional_consequence, population_freq_approx, confidence_tier, source_doi, source_detail, evidence_note
load_erap2_haplotypes <- function() {
  data <- .load_csv("erap2_haplotype_expression.csv",
                     c("haplotype", "tag_snp_rs2248374", "expression_status", "confidence_tier", "source_doi"),
                     "erap2_haplotypes")
  expected <- c("HapA", "HapB", "HapC")
  missing <- setdiff(expected, data$haplotype)
  if (length(missing) > 0) {
    stop(sprintf("[erap2_haplotypes] Missing haplotypes: %s", paste(missing, collapse = ", ")),
         call. = FALSE)
  }
  data
}

#' Load ERAP rsID -> GRCh38 VCF-allele -> amino-acid crosswalk
#'
#' Callset-agnostic (pure GRCh38) mapping that makes the ERAP1 allotype and
#' ERAP2 haplotype tables executable against a GRCh38 VCF/PGEN. ERAP1 and ERAP2
#' sit in opposite orientation at 5q15, so coding_strand is recorded per gene
#' and ref_aa/alt_aa are strand-consistent (coding strand) while ref_allele/
#' alt_allele are GRCh38 forward-strand (VCF) bases.
#' @return tibble: gene, rsid, aa_position, aa_ancestral, aa_derived, grch38_chrom,
#'   grch38_pos, ref_allele, alt_allele, ref_aa, alt_aa, coding_strand,
#'   dbsnp_build, ensembl_release, confidence_tier, source, evidence_note
load_erap_crosswalk <- function() {
  data <- .load_csv("erap_snp_crosswalk.csv",
                     c("gene", "rsid", "aa_position", "grch38_chrom", "grch38_pos",
                       "ref_allele", "alt_allele", "coding_strand", "confidence_tier", "source"),
                     "erap_crosswalk")
  bad_gene <- setdiff(unique(data$gene), c("ERAP1", "ERAP2"))
  if (length(bad_gene) > 0) {
    stop(sprintf("[erap_crosswalk] Unexpected gene values: %s", paste(bad_gene, collapse = ", ")),
         call. = FALSE)
  }
  bad_strand <- setdiff(unique(data$coding_strand), c("plus", "minus"))
  if (length(bad_strand) > 0) {
    stop(sprintf("[erap_crosswalk] coding_strand must be plus/minus, got: %s", paste(bad_strand, collapse = ", ")),
         call. = FALSE)
  }
  stopifnot("grch38_pos must be numeric" = is.numeric(data$grch38_pos))
  data
}

#' Load HLA-A lineage-level expression z-scores
#' @return tibble: hla_a_lineage, expression_z_score, expression_class, confidence_tier, source_doi, source_detail, evidence_note
load_hla_a_expression <- function() {
  data <- .load_csv("hla_a_estimated_expression.csv",
                     c("hla_a_lineage", "expression_z_score", "expression_class", "confidence_tier"),
                     "hla_a_expression")
  bad <- setdiff(unique(data$expression_class), c("high", "medium", "low"))
  if (length(bad) > 0) {
    stop(sprintf("[hla_a_expression] Unexpected expression_class: %s", paste(bad, collapse = ", ")),
         call. = FALSE)
  }
  stopifnot("expression_z_score must be numeric" = is.numeric(data$expression_z_score))
  stopifnot("All lineages must start with A*" = all(str_starts(data$hla_a_lineage, "A\\*")))
  data
}

#' Load HLA-C allotype-level expression (MFI + z-scores)
#' @return tibble: hla_c_allotype, expression_mfi, expression_z_score, expression_class, confidence_tier, source_doi, source_detail, evidence_note
load_hla_c_expression <- function() {
  data <- .load_csv("hla_c_expression.csv",
                     c("hla_c_allotype", "expression_mfi", "expression_z_score", "expression_class", "confidence_tier"),
                     "hla_c_expression")
  bad <- setdiff(unique(data$expression_class), c("high", "medium", "low"))
  if (length(bad) > 0) {
    stop(sprintf("[hla_c_expression] Unexpected expression_class: %s", paste(bad, collapse = ", ")),
         call. = FALSE)
  }
  stopifnot("expression_mfi must be numeric" = is.numeric(data$expression_mfi))
  stopifnot("expression_z_score must be numeric" = is.numeric(data$expression_z_score))
  stopifnot("All allotypes must start with C*" = all(str_starts(data$hla_c_allotype, "C\\*")))
  data
}

#' Load HLA functional divergence scores (pairwise)
#' @param locus optional; filter to "A", "B", or "C"
#' @return tibble: allele_1, allele_2, locus, fd_score, confidence_tier, source_doi
load_hla_divergence <- function(locus = NULL) {
  data <- .load_csv("hla_functional_divergence.csv",
                     c("allele_1", "allele_2", "locus", "fd_score"),
                     "hla_divergence")
  bad_locus <- setdiff(unique(data$locus), c("A", "B", "C"))
  if (length(bad_locus) > 0) {
    stop(sprintf("[hla_divergence] Unexpected locus values: %s", paste(bad_locus, collapse = ", ")),
         call. = FALSE)
  }
  stopifnot("fd_score must be numeric" = is.numeric(data$fd_score))
  if (!is.null(locus)) {
    data <- data %>% filter(locus == !!locus)
  }
  data
}

# ============================================================
# Backward-compatible vector getters
# ============================================================
# These return vector formats commonly expected by downstream
# analysis projects (classify_bw, classify_c_group, etc.)

# --- Bw4/Bw6 vectors ---
# Pipeline expects: character vectors of 2-field alleles (e.g., "B*57:01")
# Used in classify_bw() via %in%

#' Get HLA-B Bw4-80I alleles (replaces BW4_80I_ALLELES constant)
#' @return character vector of 2-field HLA-B alleles with Bw4-80I
get_bw4_80i_alleles <- function() {
  data <- load_bw4_80i_classification()
  data %>%
    filter(locus == "B", kir_ligand == "Bw4 - 80I") %>%
    pull(allele_2field)
}

#' Get HLA-B Bw4-80T alleles (replaces BW4_80T_ALLELES constant)
#' @return character vector of 2-field HLA-B alleles with Bw4-80T
get_bw4_80t_alleles <- function() {
  data <- load_bw4_80i_classification()
  data %>%
    filter(locus == "B", kir_ligand == "Bw4 - 80T") %>%
    pull(allele_2field)
}

#' Get HLA-A Bw4-80I alleles (replaces HLA_A_BW4_80I_ALLELES constant)
#' @return character vector of 2-field HLA-A alleles with Bw4-80I
get_hla_a_bw4_alleles <- function() {
  data <- load_bw4_80i_classification()
  data %>%
    filter(locus == "A", kir_ligand == "Bw4 - 80I") %>%
    pull(allele_2field)
}

# --- C1/C2 vectors ---
# Pipeline expects: character vectors of 2-field alleles (e.g., "C*07:02")
# Used in classify_c_group() via %in%

#' Get HLA-C C1 alleles (replaces C1_ALLELES constant)
#' @return character vector of 2-field HLA-C alleles with C1
get_c1_alleles <- function() {
  data <- load_c1_c2_classification()
  data %>% filter(kir_ligand == "C1") %>% pull(allele_2field)
}

#' Get HLA-C C2 alleles (replaces C2_ALLELES constant)
#' @return character vector of 2-field HLA-C alleles with C2
get_c2_alleles <- function() {
  data <- load_c1_c2_classification()
  data %>% filter(kir_ligand == "C2") %>% pull(allele_2field)
}

# --- KIR3DL1 expression vectors ---
# Pipeline expects: character vectors of 3-digit allele group codes (e.g., "001")
# Used in classify_kir3dl1() which extracts 3-digit groups from allele strings

#' Get KIR3DL1 high-expression alleles (replaces KIR3DL1_HIGH constant)
#' @return character vector of 3-digit allele codes
get_kir3dl1_high <- function() {
  data <- load_kir3dl1_expression()
  data %>% filter(expression_class == "high") %>% pull(allele)
}

#' Get KIR3DL1 low-expression alleles (replaces KIR3DL1_LOW constant)
#' @return character vector of 3-digit allele codes
get_kir3dl1_low <- function() {
  data <- load_kir3dl1_expression()
  data %>% filter(expression_class == "low") %>% pull(allele)
}

#' Get KIR3DL1 null-expression alleles (replaces KIR3DL1_NULL constant)
#' @return character vector of 3-digit allele codes
get_kir3dl1_null <- function() {
  data <- load_kir3dl1_expression()
  data %>% filter(expression_class == "null") %>% pull(allele)
}

# --- KIR3DL1 binding vectors ---
# Backward-compatible getters filter the full binding matrix to A*32:01

#' Get KIR3DL1-A*32:01 MFI named vector (replaces KIR3DL1_A3201_MFI constant)
#' @return named numeric vector: 3-digit KIR3DL1 allele code -> MFI value
get_kir3dl1_mfi <- function() {
  data <- load_kir3dl1_binding(hla = "A*32:01")
  mfi <- data$mfi_mean
  names(mfi) <- data$kir3dl1_allele
  mfi[!is.na(mfi)]
}

#' Get KIR3DL1-HLA binding matrix as a tidy tibble
#' @return tibble: kir3dl1_allele, hla_allele, mfi_mean
get_kir3dl1_binding_matrix <- function() {
  load_kir3dl1_binding() %>%
    select(kir3dl1_allele, hla_allele, mfi_mean)
}

# --- HLA-B -21 leader lookup ---
# Returns lookup from allele -> M/T (allele-level, replacing family-level assignments)
# This provides allele-level lookup which is more accurate

#' Get HLA-B -21 leader lookup vector
#' @return named character vector: 2-field allele (e.g., "B*44:02") -> "M" or "T"
get_b_leader_lookup <- function() {
  data <- load_b_leader_assignments()
  leaders <- data$b_leader
  names(leaders) <- data$allele
  leaders
}

# --- HLA-A expression lookup ---

#' Get HLA-A expression z-score named vector (replaces HLA_A_EXPRESSION_REF constant)
#' @return named numeric vector: lineage (e.g., "A*01") -> z-score
get_hla_a_expression_zscore <- function() {
  data <- load_hla_a_expression()
  zscores <- data$expression_z_score
  names(zscores) <- data$hla_a_lineage
  zscores
}

#' Get HLA-A expression classified lookup (replaces HLA_A_EXPR_CLASSIFIED constant)
#' @return named character vector: lineage (e.g., "A*01") -> "high"/"medium"/"low"
get_hla_a_expression_class <- function() {
  data <- load_hla_a_expression()
  classes <- data$expression_class
  names(classes) <- data$hla_a_lineage
  classes
}

# --- HLA-C expression lookup ---

#' Get HLA-C expression MFI named vector
#' @return named numeric vector: allotype (e.g., "C*07") -> MFI
get_hla_c_expression_mfi <- function() {
  data <- load_hla_c_expression()
  mfi <- data$expression_mfi
  names(mfi) <- data$hla_c_allotype
  mfi
}

#' Get HLA-C expression classified lookup
#' @return named character vector: allotype (e.g., "C*07") -> "high"/"medium"/"low"
get_hla_c_expression_class <- function() {
  data <- load_hla_c_expression()
  classes <- data$expression_class
  names(classes) <- data$hla_c_allotype
  classes
}

# --- ERAP crosswalk lookup ---

#' Get ERAP SNP crosswalk as a GRCh38 VCF-coordinate lookup
#' @return tibble: rsid, gene, grch38_chrom, grch38_pos, ref_allele, alt_allele,
#'   aa_position, ref_aa, alt_aa, coding_strand
get_erap_crosswalk_lookup <- function() {
  load_erap_crosswalk() %>%
    select(rsid, gene, grch38_chrom, grch38_pos, ref_allele, alt_allele,
           aa_position, ref_aa, alt_aa, coding_strand)
}

# --- HLA functional divergence lookup ---

#' Get HLA-A functional divergence lookup in compact pair format
#' @return tibble with geno (e.g., "A0101_A0201") and fd columns
get_hla_a_divergence_lookup <- function() {
  data <- load_hla_divergence(locus = "A")
  data %>%
    mutate(geno = paste0(
      str_replace_all(allele_1, "[*:]", ""), "_",
      str_replace_all(allele_2, "[*:]", "")
    )) %>%
    transmute(geno, fd = fd_score) %>%
    filter(!is.na(fd))
}

# ============================================================
# Legacy aliases (for existing code that uses old names)
# ============================================================

load_bw4_bw6        <- load_bw4_classification
load_bw4_80i        <- load_bw4_80i_classification
load_c1_c2          <- load_c1_c2_classification
load_b_leader       <- load_b_leader_assignments
load_erap1_activity <- load_erap1_allotypes
load_erap2_expression <- load_erap2_haplotypes
load_hla_a_expr      <- load_hla_a_expression
load_hla_c_expr      <- load_hla_c_expression

extract_bw4_vectors <- function() {
  list(
    bw4_80i      = get_bw4_80i_alleles(),
    bw4_80t      = get_bw4_80t_alleles(),
    hla_a_bw4_80i = get_hla_a_bw4_alleles()
  )
}

extract_c1c2_vectors <- function() {
  list(c1 = get_c1_alleles(), c2 = get_c2_alleles())
}

extract_kir3dl1_expression_vectors <- function() {
  list(high = get_kir3dl1_high(), low = get_kir3dl1_low(), null = get_kir3dl1_null())
}

extract_kir3dl1_mfi <- get_kir3dl1_mfi
extract_kir3dl1_binding_matrix <- get_kir3dl1_binding_matrix
