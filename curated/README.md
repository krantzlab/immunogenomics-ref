# Curated Reference Data

Literature-sourced biological classifications maintained manually. Each file includes provenance columns (`confidence_tier`, `source_doi`, `source_detail`, `evidence_note`) for full traceability. See also `../provenance.yaml`.

## Confidence Tiers

| Tier | Meaning |
|------|---------|
| A | Direct measurement for this specific allele/allotype in the cited study |
| B | Literature consensus from multiple studies or well-established classification |
| C | Inferred from subfamily similarity, estimated, or limited direct data |
| D | Single unreplicated study or expert opinion only |

## Files

| File | Rows | Content | Primary DOI |
|------|------|---------|-------------|
| `kir3dl1_expression.csv` | 17 | KIR3DL1 allele expression (high/low/null) | 10.4049/jimmunol.180.10.6743 |
| `kir3dl1_hla_binding.csv` | 323 | KIR3DL1 allotype binding to HLA class I (full MFI matrix) | 10.1101/2024.05.03.592082 |
| `erap1_allotype_activity.csv` | 10 | ERAP1 allotype trimming activity with SNP haplotypes | 10.1016/j.jbc.2021.100443 |
| `erap2_haplotype_expression.csv` | 3 | ERAP2 haplotype expression (HapA/HapB/HapC) | 10.1371/journal.pgen.1001029 |
| `hla_functional_divergence.csv` | 3,002 | HLA-A/B/C pairwise functional divergence scores | 10.1126/science.ado8609 |
| `hla_a_estimated_expression.csv` | 21 | HLA-A lineage expression z-scores (high/medium/low) | 10.1126/science.aam8825 |
| `hla_c_expression.csv` | 14 | HLA-C allotype expression MFI (high/medium/low) | 10.1182/blood-2014-09-599969 |

## Review process

1. When new literature updates a classification, edit the CSV here
2. Document the change in `provenance.yaml` (update `last_reviewed`, add reference)
3. Copy updated file to `datasets/`: `cp curated/*.csv datasets/`
4. Run validation: `Rscript -e 'source("R/validate_reference_data.R"); validate_all()'`
5. Update `CHANGELOG.md`

## Key classification decisions

These are documented in `provenance.yaml` under `key_decisions`:

- **KIR3DL1*005**: HIGH expression (ER-retained but high surface; Boudreau 2016), STRONG binding
- **KIR3DL1*004**: LOW expression (not null), WEAK binding
- **KIR3DL1*006/*007**: MFI estimated from subfamily (not directly measured); tier C
- **ERAP1 allotype 10**: very_low activity (up to 60-fold lower; protective for AS/Behcets)
- **ERAP2 HapB**: non-expressing (aberrant splicing -> NMD; maintained by balancing selection)
- **HLA-A Bw4**: Only A*23:01, A*24:02, A*32:01 are functional KIR3DL1 ligands (A*25:01 excluded per Foley 2008)
