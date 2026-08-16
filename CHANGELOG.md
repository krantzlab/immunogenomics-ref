# Changelog

## v1.2.0 (2026-08-16)

Additive release. No existing column, row or value changes, so a consumer pinned to v1.1.0 can move to v1.2.0 without altering how it reads these files.

### Managed data (database-derived)

- **`ipd_version` added to `bw4_bw6_classification.csv` and `c1_c2_classification.csv`** — the two API-derived tables recorded only `source` and `fetch_date`, so neither could state which IPD-IMGT/HLA release produced its classifications. `b_leader_assignments.csv` and `bw4_80i_classification.csv` already carried the column; all four managed tables now agree on `IPD-IMGT/HLA 3.63.0`. Column order matches the existing convention (`... source, ipd_version, fetch_date`).

  Backfilled from the recorded provenance rather than re-fetched: re-running the fetch would pull whatever IPD has released since February and change classifications, which is a separate change. Every pre-existing value is byte-identical — removing the new column reproduces the previous files exactly.

  **Additive change.** No existing column, row or value is altered, so consumers reading these files by column name are unaffected.

### Fixed

- **`CITATION.cff` declared version 1.0.0 at the v1.1.0 tag**, so anyone citing the dataset from that file cited a version that did not describe what they had. The v1.1.0 tag necessarily keeps the stale file — moving a published tag would invalidate the provenance record of every consumer pinned to it — so the correction takes effect here. The abstract also now mentions the ERAP SNP crosswalk, which v1.1.0 introduced.
- **Documented check count was stale.** `validate_all()` has registered 13 checks since `erap_crosswalk` was added in v1.1.0, but `CLAUDE.md`, `CONTRIBUTING.md` and the README's example output all still said 12. The README's curated dataset table and PASS listing were missing the crosswalk entirely.

### Infrastructure

- **Continuous integration** — `.github/workflows/validate.yml` runs `validate_all()` on every push, pull request and tag. A second, tag-only job asserts that `CITATION.cff` and `CHANGELOG.md` agree with the tag being pushed. Before this, the validation documented as a pre-commit requirement was an honour system.
- **AI attribution policy enforced mechanically.** The policy has been documented since v1.0.0 and this repository's history is clean, but that was luck of which file an assistant read — the identical policy was absent from `krantzlab/bridgie`, where 15 commits acquired `Co-Authored-By: Claude` and `Claude-Session:` trailers. New `.githooks/commit-msg` rejects them locally (`git config core.hooksPath .githooks`, once per clone); the `attribution` CI job scans every commit on a pull request for the same patterns, where neither opt-out nor `--no-verify` applies. Human `Co-Authored-By` lines are unaffected.
- **`tools/validate.R`** — single entry point for validation, used by both CI and `pixi run validate`. It exists because `validate_all()` returns its result invisibly and never signals a condition, so calling it from a shell always exited 0 even when checks failed.
- **`pixi.toml` / `pixi.lock`** — declared development environment. `default` carries what loading and validating require; `managed` adds the IPD fetch toolchain. `r-base` is pinned at 4.5.\*. HLAtools has no conda-forge build and installs into a project-local `.Rlib` via `pixi run -e managed deps`.
- `managed/fetch_kir_ligand.R` — writes `ipd_version` into both outputs. The value was already resolved from the HLAtools alignment and printed to stdout, but was never carried into the tables. An unresolved version is now written as `NA` rather than silently dropping the column.
- `R/validate_reference_data.R` — added `.check_ipd_version()`, applied to all four managed tables: the column must be present, non-blank on every row, and hold exactly one version per table. Check count is unchanged at 13; the assertion runs inside the existing per-table validators.
- `R/load_reference_data.R` — `@return` documentation for `load_bw4_classification()` and `load_c1_c2_classification()` updated.
- `provenance.yaml` — `columns` updated for both datasets.

### Documentation

- **Licence scoped to what this repository actually holds.** `LICENSE` declared MIT over the whole repository, which claimed more than the project can grant: the managed tables are derived from IPD-IMGT/HLA, and four curated tables are extracted from publications — 3,002 rows from a *Science* supplementary table among them — whose copyright is held by the publishers. MIT now covers the code and the compilation; new `DATA_TERMS.md` gives the per-file breakdown, the IPD citation and licence link, and guidance for anyone redistributing or vendoring a snapshot. No data changed and nothing became more restrictive; the previous statement was simply inaccurate.
- **`LICENSE` copyright holder filled in.** It read `Copyright (c) 2026 [Your Name]` — the template placeholder had never been replaced, leaving the grantor of the licence unstated on a public repository.
- **Release checklist** in `CONTRIBUTING.md`, recording the constraint that makes tags immutable here: downstream packages pin them and vendor `datasets/`, so a published tag is part of their provenance record and can only be superseded, never moved. It also states which schema changes are breaking.
- `build:` and `ci:` added to the commit type list, which had no category for tooling.
- README documents the pixi workflow and carries a CI status badge.

## v1.1.0 (2026-06-09)

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
