# derive_bw4_80i.R -- Derive Bw4/Bw6 and position 80 classification from protein sequences
#
# Produces:
#   datasets/bw4_80i_classification.csv
#
# Source: HLAtools protein alignment of HLA-B and HLA-A from IPD-IMGT/HLA
#
# Method:
#   - Position 83: R = Bw4, G = Bw6 (primary Bw4/Bw6 discriminator)
#                  anything else = Non-canonical: the allele carries neither
#                  canonical epitope. R83 is the residue KIR3DL1 binding
#                  requires, so this is a determinate finding, not missing data.
#   - Position 80: I (isoleucine) = strong KIR3DL1 ligand
#                  T (threonine)  = weaker KIR3DL1 ligand
#                  N (asparagine) = Bw6 (not a KIR3DL1 ligand)
#
# Null alleles (N suffix) are excluded from both loci. They are not expressed,
# so they present no epitope and cannot be KIR ligands; and since they are the
# only source of unreadable sequence, excluding them leaves every row with a
# readable position 83. That is what lets "not a Bw4/Bw6 call" carry exactly
# one meaning. See classification_status.
#
# Includes HLA-A alleles that carry the Bw4 epitope and function as
# KIR3DL1 ligands (A*23:01, A*24:02, A*32:01).
# Note: A*25:01 has Bw4 serologically but does NOT inhibit KIR3DL1+ NK cells
# (Foley et al., Blood 2008; DOI: 10.1182/blood-2008-03-142943).

library(HLAtools)
library(dplyr)
library(stringr)
library(readr)
library(here)

cat("=== derive_bw4_80i.R ===\n")

# ============================================================
# 1. Build protein alignments from IMGT/HLA via HLAtools
# ============================================================
cat("  Building protein alignments from IMGT/HLA ...\n")

ipd_version <- NULL
B_align <- NULL
A_align <- NULL

for (v in c("3.63.0", "3.62.0", "3.61.0", "3.60.0", "3.58.0", "3.56.0")) {
  B_align <- tryCatch(buildAlignments("B", source = "AA", version = v),
                       error = function(e) NULL)
  if (!is.null(B_align)) {
    ipd_version <- B_align$B$Version
    A_align <- tryCatch(buildAlignments("A", source = "AA", version = v),
                         error = function(e) NULL)
    cat(sprintf("  Alignment version: %s\n", ipd_version))
    break
  }
}

if (is.null(B_align)) stop("Could not build HLA-B protein alignment from HLAtools")

# ============================================================
# 2. Classify HLA-B alleles
# ============================================================
cat("  Classifying HLA-B alleles ...\n")

B_aa <- B_align$B$AA

hla_b_classified <- B_aa %>%
  select(allele, trimmed_allele, pos77 = `77`, pos78 = `78`, pos79 = `79`,
         pos80 = `80`, pos81 = `81`, pos82 = `82`, pos83 = `83`) %>%
  filter(!duplicated(trimmed_allele)) %>%
  # Null alleles are excluded. A null allele is not expressed, presents no
  # epitope, and therefore cannot be a KIR ligand -- and reading positions
  # 77-83 out of a frameshifted or absent sequence is what produced the
  # spurious Bw4-80D and Bw4-80E categories in v1.x (B*35:542N read as
  # "EPAEPAR", which ends in R and so passed the position-83 rule).
  #
  # It also keeps the table's meaning single-valued: every remaining row has a
  # readable sequence, so "not classified" can only mean "position 83 is
  # neither R nor G", never "no sequence to read".
  filter(!str_detect(trimmed_allele, "N$")) %>%
  mutate(
    bw_group = case_when(
      pos83 == "R" ~ "Bw4",
      pos83 == "G" ~ "Bw6",
      TRUE ~ "non_canonical"
    ),
    pos80_class = case_when(
      bw_group == "Bw4" & pos80 == "I" ~ "80I",
      bw_group == "Bw4" & pos80 == "T" ~ "80T",
      bw_group == "Bw4" ~ paste0("80", pos80),  # rare non-I/T variants
      TRUE ~ NA_character_
    ),
    kir_ligand = case_when(
      bw_group == "Bw6" ~ "Bw6",
      bw_group == "Bw4" & pos80 == "I" ~ "Bw4 - 80I",
      bw_group == "Bw4" & pos80 == "T" ~ "Bw4 - 80T",
      bw_group == "Bw4" ~ paste0("Bw4 - ", pos80_class),
      TRUE ~ "Non-canonical"
    ),
    # kir_ligand is the display value and keeps its historical spacing.
    # kir_ligand_code is the identifier: it is what ends up in a feature id, a
    # column name or a model matrix, so it carries no spaces or punctuation.
    kir_ligand_code = case_when(
      kir_ligand == "Bw6" ~ "Bw6",
      kir_ligand == "Non-canonical" ~ "NonCanonical",
      TRUE ~ paste0("Bw4_", pos80_class)
    ),
    # Records why a row is not a Bw4/Bw6 call. Only one reason is now possible,
    # because null alleles -- the sole source of unreadable sequence -- are
    # filtered above.
    classification_status = if_else(bw_group == "non_canonical",
                                    "non_canonical_83", "classified"),
    bw4_motif = paste0(pos77, pos78, pos79, pos80, pos81, pos82, pos83)
  ) %>%
  transmute(
    allele_2field = trimmed_allele,
    locus = "B",
    kir_ligand,
    kir_ligand_code,
    classification_status,
    pos80 = pos80,
    bw4_motif_77_83 = bw4_motif,
    source = "HLAtools_protein_alignment",
    ipd_version = ipd_version,
    fetch_date = as.character(Sys.time())
  )

cat(sprintf("  HLA-B: %d unique 2-field alleles classified\n", nrow(hla_b_classified)))
cat(sprintf("  Distribution:\n"))
print(table(hla_b_classified$kir_ligand))

# ============================================================
# 3. Classify HLA-A alleles (Bw4 carriers)
# ============================================================
cat("\n  Classifying HLA-A alleles ...\n")

hla_a_classified <- tibble()

if (!is.null(A_align)) {
  A_aa <- A_align$A$AA

  hla_a_raw <- A_aa %>%
    select(allele, trimmed_allele, pos80 = `80`, pos83 = `83`) %>%
    filter(!duplicated(trimmed_allele))

  # HLA-A alleles with Bw4: check position 83 for R (same as HLA-B)
  # Also check known functional KIR3DL1 ligands
  hla_a_classified <- hla_a_raw %>%
    filter(pos83 == "R" | trimmed_allele %in% c("A*23:01", "A*24:02", "A*32:01")) %>%
    # Null alleles excluded for the same reason as at HLA-B.
    filter(!str_detect(trimmed_allele, "N$")) %>%
    mutate(
      kir_ligand = case_when(
        pos80 == "I" ~ "Bw4 - 80I",
        pos80 == "T" ~ "Bw4 - 80T",
        TRUE ~ paste0("Bw4 - 80", pos80)
      ),
      kir_ligand_code = paste0("Bw4_80", pos80),
      # Every HLA-A row reaches this table by carrying Bw4, so there is no
      # non-canonical case here: the filter above selects on position 83.
      classification_status = "classified"
    ) %>%
    transmute(
      allele_2field = trimmed_allele,
      locus = "A",
      kir_ligand,
      kir_ligand_code,
      classification_status,
      pos80,
      bw4_motif_77_83 = NA_character_,
      source = "HLAtools_protein_alignment",
      ipd_version = ipd_version,
      fetch_date = as.character(Sys.time())
    )

  cat(sprintf("  HLA-A Bw4 alleles found: %d\n", nrow(hla_a_classified)))
  if (nrow(hla_a_classified) > 0) {
    cat("  HLA-A Bw4 alleles:\n")
    hla_a_classified %>%
      select(allele_2field, kir_ligand, pos80) %>%
      as.data.frame() %>%
      print(row.names = FALSE)
  }

  # Note about A*25:01
  a2501 <- hla_a_raw %>% filter(trimmed_allele == "A*25:01")
  if (nrow(a2501) > 0) {
    cat(sprintf("\n  Note: A*25:01 pos80=%s, pos83=%s - serological Bw4 but\n",
                a2501$pos80[1], a2501$pos83[1]))
    cat("    does NOT inhibit KIR3DL1+ NK cells (Foley 2008); included in output\n")
    cat("    but downstream analyses should handle exclusion based on functional data.\n")
  }
}

# ============================================================
# 4. Combine and save
# ============================================================
bw4_combined <- bind_rows(hla_b_classified, hla_a_classified) %>%
  arrange(locus, allele_2field)

# ------------------------------------------------------------
# 4a. Carry the API's own call alongside the sequence-derived one
# ------------------------------------------------------------
# The IPD Allele Query API classifies some alleles that this alignment calls
# non-canonical -- 24 of them at 3.63.0, all expressed. Which source is right
# is a scientific question about whether IPD's curated kir_ligand field or a
# strict position-83 rule better reflects KIR3DL1 biology, and it is not
# settled here. Recording both makes the disagreement visible and countable
# instead of a silent editorial choice; a jump in the count across IPD
# releases means something changed upstream.
#
# HLA-A is absent from the API table by design, so api_kir_ligand is NA there.
api_path <- here("datasets", "bw4_bw6_classification.csv")
if (file.exists(api_path)) {
  api_calls <- read_csv(api_path, show_col_types = FALSE) %>%
    select(allele_2field, api_kir_ligand = kir_ligand)
  bw4_combined <- bw4_combined %>%
    left_join(api_calls, by = "allele_2field") %>%
    relocate(api_kir_ligand, .after = bw4_motif_77_83)

  n_dis <- sum(!is.na(bw4_combined$api_kir_ligand) &
                 bw4_combined$api_kir_ligand != bw4_combined$kir_ligand)
  cat(sprintf("\n  API/sequence disagreements recorded: %d\n", n_dis))
} else {
  warning("bw4_bw6_classification.csv not found; api_kir_ligand will be absent. ",
          "Run fetch_kir_ligand.R first.", call. = FALSE)
  bw4_combined$api_kir_ligand <- NA_character_
  bw4_combined <- bw4_combined %>% relocate(api_kir_ligand, .after = bw4_motif_77_83)
}

# `allele` is the uniform key across every managed table; `allele_2field` is
# kept as an exact duplicate for consumers written against the old name, and is
# deprecated for removal in v3.0.0.
bw4_combined$allele <- bw4_combined$allele_2field
bw4_combined <- bw4_combined %>% relocate(allele, .before = allele_2field)

output_path <- here("datasets", "bw4_80i_classification.csv")
write_csv(bw4_combined, output_path)

cat(sprintf("\n  Wrote %s (%d rows)\n", basename(output_path), nrow(bw4_combined)))
cat(sprintf("    HLA-B: %d (Bw4-80I=%d, Bw4-80T=%d, Bw6=%d)\n",
            sum(bw4_combined$locus == "B"),
            sum(bw4_combined$locus == "B" & bw4_combined$kir_ligand == "Bw4 - 80I"),
            sum(bw4_combined$locus == "B" & bw4_combined$kir_ligand == "Bw4 - 80T"),
            sum(bw4_combined$locus == "B" & bw4_combined$kir_ligand == "Bw6")))
cat(sprintf("    HLA-A: %d (all Bw4)\n", sum(bw4_combined$locus == "A")))

# ============================================================
# 5. Compare with previous hardcoded classifications
# ============================================================
cat("\n  Comparison with previous hardcoded classifications:\n")

old_80i <- c("B*27:02", "B*27:05", "B*38:01", "B*44:02", "B*44:03",
             "B*47:01", "B*49:01", "B*51:01", "B*51:02", "B*52:01",
             "B*53:01", "B*57:01", "B*57:03", "B*58:01", "B*58:02", "B*59:01")

for (a in old_80i) {
  row <- bw4_combined %>% filter(allele_2field == a)
  if (nrow(row) > 0) {
    new_class <- row$kir_ligand[1]
    status <- ifelse(new_class == "Bw4 - 80I", "OK", sprintf("CHANGED → %s", new_class))
    cat(sprintf("    %s: %s\n", a, status))
  }
}

cat(sprintf("\n  IPD-IMGT/HLA version: %s\n", ipd_version))
cat(sprintf("  Derivation date: %s\n", as.character(Sys.time())))
