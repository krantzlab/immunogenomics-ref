# Managed Reference Data

Scripts that query authoritative databases to generate reference lookup tables.

## Scripts

| Script | Primary Source | Output |
|--------|---------------|--------|
| `fetch_kir_ligand.R` | IPD Allele Query API + HLAtools cross-validation | `datasets/bw4_bw6_classification.csv`, `datasets/c1_c2_classification.csv` |
| `fetch_b_leader.R` | HLAtools protein alignment (position -21) + API cross-validation | `datasets/b_leader_assignments.csv` |
| `derive_bw4_80i.R` | HLAtools protein alignment (positions 77-83, 80) | `datasets/bw4_80i_classification.csv` |

## How to update

Run scripts in order:

```bash
Rscript managed/fetch_kir_ligand.R
Rscript managed/fetch_b_leader.R
Rscript managed/derive_bw4_80i.R
```

Requirements: `httr`, `jsonlite`, `HLAtools`, `dplyr`, `stringr`, `readr`, `here`.

All scripts automatically:
- Cross-validate between API and sequence-derived classifications
- Log IPD-IMGT/HLA version and fetch date
- Report any discrepancies

## Data sources

- **IPD Allele Query API**: `https://www.ebi.ac.uk/cgi-bin/ipd/api/allele` (paginated, max 1000/page)
- **HLAtools**: R package wrapping IMGT/HLA protein alignments from ANHIG/IMGTHLA GitHub

## Last fetch

- **Date**: 2026-02-19
- **IPD-IMGT/HLA version**: 3.63.0
- **Results**:
  - `bw4_bw6_classification.csv`: 6,350 HLA-B alleles (Bw4-80I: 1,231, Bw4-80T: 1,157, Bw6: 3,957)
  - `c1_c2_classification.csv`: 4,978 HLA-C alleles (C1: 3,114, C2: 1,864)
  - `b_leader_assignments.csv`: 4,353 HLA-B alleles (M: 1,006, T: 3,347)
  - `bw4_80i_classification.csv`: 7,953 alleles (6,848 HLA-B + 1,105 HLA-A)
- **Cross-validation**: 6,258 API/sequence matches, 11 mismatches corrected; 0 C1/C2 mismatches; 0 leader mismatches

## Key findings from sequence data

4 alleles reclassified from previous hardcoded Bw4-80I to Bw4-80T:
- B\*27:05, B\*44:02, B\*44:03, B\*47:01

6 HLA-B families reclassified from Petersdorf 2022 M-leader to T-leader:
- B\*44, B\*45, B\*47, B\*56, B\*57, B\*58 (only B\*07 remains M)
