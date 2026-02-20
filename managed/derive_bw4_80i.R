# derive_bw4_80i.R -- Derive Bw4/Bw6 and position 80 classification from protein sequences
#
# Produces:
#   datasets/bw4_80i_classification.csv
#
# Source: HLAtools protein alignment of HLA-B and HLA-A from IPD-IMGT/HLA
#
# Method:
#   - Position 83: R = Bw4, G = Bw6 (primary Bw4/Bw6 discriminator)
#   - Position 80: I (isoleucine) = strong KIR3DL1 ligand
#                  T (threonine)  = weaker KIR3DL1 ligand
#                  N (asparagine) = Bw6 (not a KIR3DL1 ligand)
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
  mutate(
    bw_group = case_when(
      pos83 == "R" ~ "Bw4",
      pos83 == "G" ~ "Bw6",
      TRUE ~ "unclassified"
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
      TRUE ~ "unclassified"
    ),
    bw4_motif = paste0(pos77, pos78, pos79, pos80, pos81, pos82, pos83)
  ) %>%
  transmute(
    allele_2field = trimmed_allele,
    locus = "B",
    kir_ligand,
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
    mutate(
      kir_ligand = case_when(
        pos80 == "I" ~ "Bw4 - 80I",
        pos80 == "T" ~ "Bw4 - 80T",
        TRUE ~ paste0("Bw4 - 80", pos80)
      )
    ) %>%
    transmute(
      allele_2field = trimmed_allele,
      locus = "A",
      kir_ligand,
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
