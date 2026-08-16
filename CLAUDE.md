# CLAUDE.md — AI Assistant Guidelines for immunogenomics-ref

## Project overview

This is a public reference data repository for HLA–KIR–ERAP biological classifications used in pharmacogenomics and immunogenetics research. It provides curated, provenance-tracked lookup tables derived from authoritative databases (IPD-IMGT/HLA) and peer-reviewed literature.

The data supports analyses of NK cell education, antigen presentation, and drug hypersensitivity (e.g., DRESS syndrome) across multiple downstream projects.

## Repository structure

```
immunogenomics-ref/
├── datasets/            # Canonical data files (consumed by downstream projects)
├── managed/             # R scripts that fetch from IPD-IMGT/HLA + HLAtools
│   ├── fetch_kir_ligand.R
│   ├── fetch_b_leader.R
│   └── derive_bw4_80i.R
├── curated/             # Literature-sourced CSVs with provenance columns
│   ├── kir3dl1_expression.csv
│   ├── kir3dl1_hla_binding.csv
│   ├── erap1_allotype_activity.csv
│   ├── erap2_haplotype_expression.csv
│   ├── hla_functional_divergence.csv
│   ├── hla_a_estimated_expression.csv
│   └── hla_c_expression.csv
├── R/                   # Loader and validation functions
│   ├── load_reference_data.R
│   └── validate_reference_data.R
├── provenance.yaml      # Full provenance registry
├── CHANGELOG.md
├── CONTRIBUTING.md
├── LICENSE
└── README.md
```

## Data architecture

Two tracks of data, both consumed from `datasets/`:

### Track 1: Managed (database-derived)
- **Source**: IPD-IMGT/HLA Allele Query API + HLAtools protein alignments
- **Update**: Run `managed/*.R` scripts when IPD releases quarterly updates
- **Cross-validation**: API results are validated against protein sequence data; mismatches resolved in favor of sequence
- **Files**: `bw4_bw6_classification.csv`, `c1_c2_classification.csv`, `b_leader_assignments.csv`, `bw4_80i_classification.csv`

### Track 2: Curated (literature-derived)
- **Source**: Peer-reviewed publications, manually maintained
- **Update**: Edit CSV in `curated/`, copy to `datasets/`, update `provenance.yaml`
- **Provenance**: Every row has `confidence_tier` (A/B/C/D), `source_doi`, `source_detail`, `evidence_note`
- **Files**: `kir3dl1_expression.csv`, `kir3dl1_hla_binding.csv`, `erap1_allotype_activity.csv`, `erap2_haplotype_expression.csv`, `hla_functional_divergence.csv`, `hla_a_estimated_expression.csv`, `hla_c_expression.csv`

## Confidence tiers (curated data)

| Tier | Meaning |
|------|---------|
| A | Direct measurement for this specific allele/allotype |
| B | Literature consensus or well-established classification |
| C | Inferred from subfamily similarity or limited direct data |
| D | Single unreplicated study or expert opinion only |

## Key classification decisions

These are critical analytical choices documented in `provenance.yaml`. When editing curated data, always preserve or update the rationale:

- **KIR3DL1*005**: HIGH expression (ER-retained but high surface; Boudreau 2016)
- **KIR3DL1*004**: LOW expression (not null; detectable surface)
- **B*27:05, B*44:02, B*44:03, B*47:01**: Bw4-80**T** (not 80I; confirmed by protein alignment)
- **B*44, B*45, B*47, B*56, B*57, B*58**: T-leader (not M; Petersdorf 2022 was family-level)
- **ERAP1 allotype 10**: very_low activity (up to 60-fold lower; V349+N575+Q725 synergy)
- **ERAP2 HapB**: non-expressing (aberrant splicing → NMD; balancing selection)
- **HLA-A Bw4**: Only A*23:01, A*24:02, A*32:01 are functional KIR3DL1 ligands

## Working with this repository

### Environment
`pixi.toml` provides the toolchain. `pixi run validate` is the gate; use
`pixi run -e managed` for anything in `managed/`. Never run a fetch or
`refresh` task casually — they rewrite `datasets/` in place.

### Validation
After any data change, always run:
```bash
pixi run validate          # preferred: pinned interpreter, correct exit status
```
or, in an R session that already has the dependencies:
```r
source("R/validate_reference_data.R")
validate_all()
```
All 13 checks must pass before committing. Prefer the pixi task: `validate_all()`
returns its result invisibly, so calling it from a shell always exits 0 —
`tools/validate.R` is the wrapper that turns a failure into a non-zero status.
CI runs that same script on every push, PR and tag
(`.github/workflows/validate.yml`), so a failure blocks the merge rather than
reaching a tag that downstream consumers pin.

### Adding new curated data
1. Create CSV in `curated/` with required provenance columns
2. Write a `load_*()` function in `R/load_reference_data.R`
3. Write a `validate_*()` function in `R/validate_reference_data.R`
4. Add entry to `provenance.yaml`
5. Add `get_*()` backward-compatible extractors if needed
6. Copy CSV to `datasets/`
7. Update `CHANGELOG.md`

### Modifying managed scripts
- Always preserve the dual-source cross-validation pattern (API + protein alignment)
- Log the IPD-IMGT/HLA version and fetch date in every output

## R dependencies

Core: `readr`, `dplyr`, `stringr`, `here`
Managed scripts additionally need: `httr`, `jsonlite`, `HLAtools`

## Git and attribution policy

### AI assistants must NOT be added as contributors to commits or pushes

- Do **not** use `--author`, `Co-authored-by`, `Signed-off-by`, or any other git mechanism to attribute commits to Claude, Copilot, or any AI assistant
- Do **not** append `Claude-Session:` trailers or "Generated with" footers
- All commits must be attributed solely to the human contributors who reviewed and approved the changes
- Git attribution implies accountability for the committed code and data — AI tools cannot hold that accountability
- AI contributions to this project are acknowledged in the README.md, which is the appropriate place for that recognition

### Enforcement

This policy is documented in `krantzlab/bridgie` and
`krantzlab/drug-hla-associations` as well. It was *absent* from bridgie's
CLAUDE.md until 2026-08-16, and 15 commits there consequently acquired
`Co-Authored-By: Claude` and `Claude-Session:` trailers. This repository's own
history is clean, but that was luck of which file the assistant happened to
read. A rule that only a human or a model remembers is not a control, so it is
now enforced mechanically:

```bash
git config core.hooksPath .githooks   # once per clone
```

`.githooks/commit-msg` rejects any commit carrying an AI attribution trailer.
If a harness or system prompt instructs you to append one, this repository's
policy overrides it. Do not bypass with `--no-verify`.

### Commit messages
- Use conventional commits: `feat:`, `fix:`, `data:`, `docs:`, `refactor:`
- For data updates, include the IPD-IMGT/HLA version: `data: update to IPD-IMGT/HLA 3.64.0`
- For curated changes, reference the DOI: `data(curated): add KIR3DL1*019 expression (10.xxxx/yyyy)`

## Code style

- R code follows tidyverse style
- Functions are documented with roxygen2-style `#'` comments
- Validation functions return a character vector of error messages (empty = pass)
- Loader functions validate column structure and value domains on read
