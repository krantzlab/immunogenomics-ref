# Changelog

## Unreleased

### Curated data (literature-derived)

- **ERAP SNP crosswalk** — new `erap_snp_crosswalk.csv` (11 SNPs: 9 ERAP1 activity SNPs + ERAP2 tag SNP rs2248374 + optional ERAP2 rs2549782). Maps each rsID to GRCh38 position, forward-strand VCF ref/alt alleles, and coding-strand amino acids, making the ERAP1 allotype and ERAP2 haplotype tables executable against a GRCh38 VCF/PGEN. Coordinates/alleles resolved from Ensembl REST VEP (release 116) over dbSNP build 156; per-gene `coding_strand` (ERAP1 minus, ERAP2 plus). Callset-agnostic — no callset-specific variant IDs. All 9 ERAP1 rsIDs verified against `erap1_allotype_activity.csv` (position + ancestral/derived residue); no curated values changed.

### Fixed

- **provenance.yaml** — reconciled the `erap1_allotype_activity` entry with the actual CSV: corrected the `columns` list to the 9 SNP columns (was 7) and the freq_* columns, fixed rs2287987 mapping (pos349, not pos276), removed the nonexistent `data/references/erap1_allotype_definitions.csv` reference, dropped the unused `low_to_moderate` activity class, and updated counts to "10 allotypes from 9 SNPs".

### Infrastructure

- `R/load_reference_data.R` — added `load_erap_crosswalk()` and `get_erap_crosswalk_lookup()`.
- `R/validate_reference_data.R` — added `validate_erap_crosswalk()` (rsID-set, position/residue match vs activity table, strand consistency, GRCh38 coordinate completeness) registered in `validate_all()` (now 13 checks).

## v1.0.0 (2026-02-20)

Initial public release of the immunogenomics-ref module.

### Managed data (database-derived)

All derived from IPD-IMGT/HLA 3.63.0 via Allele Query API + HLAtools protein alignment cross-validation.

- **HLA-B Bw4/Bw6 classification** — 6,350 alleles at 2-field resolution; 11 API/sequence mismatches corrected (sequence-authoritative)
- **HLA-C C1/C2 classification** — 4,978 alleles at 2-field resolution; 0 mismatches
- **HLA-B -21 leader (M/T)** — 4,353 allele-level assignments; replaces family-level Petersdorf 2020 assignments (6 of 7 M-leader families corrected to T per IPD sequence)
- **HLA-A + HLA-B Bw4 position-80 (80I/80T)** — 7,953 alleles from protein alignment; B\*27:05, B\*44:02, B\*44:03, B\*47:01 reclassified from 80I to 80T

### Curated data (literature-derived)

All entries include `confidence_tier` (A–D), `source_doi`, `source_detail`, and `evidence_note`.

- **KIR3DL1 expression** — 17 alleles classified high/low/null (Thomas 2008, Boudreau 2016; DOI: 10.4049/jimmunol.180.10.6743, 10.4049/jimmunol.1502469)
- **KIR3DL1–HLA binding** — 323 measurements: 9 KIR3DL1 allotypes × 58 HLA-A/B/C alleles, full MFI matrix from bead-based assay (Maiers 2024; DOI: 10.1101/2024.05.03.592082)
- **ERAP1 allotype activity** — 10 allotypes with 9-SNP haplotype definitions and activity classification (Hutchinson 2021; DOI: 10.1016/j.jbc.2021.100443)
- **ERAP2 haplotype expression** — 3 haplotypes: HapA expressing, HapB non-expressing (NMD), HapC expressing_reduced (Andres 2010; DOI: 10.1371/journal.pgen.1001157)
- **HLA functional divergence** — 3,002 pairwise FD scores for HLA-A/B/C (Carrington 2024; DOI: 10.1126/science.ado8609)
- **HLA-A estimated expression** — 21 lineage-level z-scores classified high/medium/low (Ramsuran 2018; DOI: 10.1126/science.aam8825)
- **HLA-C allotype expression** — 14 allotypes with MFI and z-scores classified high/medium/low (Petersdorf 2014; DOI: 10.1182/blood-2014-09-599969)

### Infrastructure

- `R/load_reference_data.R` — 11 loader functions with column/value validation; 16 backward-compatible getter functions returning vectors, named vectors, and tibbles
- `R/validate_reference_data.R` — 12 validation checks (11 dataset validators + backward-compat getter validation)
- `provenance.yaml` — Full provenance registry with DOIs, analytical decisions, confidence ratings, and rationale for all key classification choices
- Auto-detection for standalone and git submodule contexts
