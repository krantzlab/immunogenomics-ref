# immunogenomics-ref

[![validate](https://github.com/krantzlab/immunogenomics-ref/actions/workflows/validate.yml/badge.svg)](https://github.com/krantzlab/immunogenomics-ref/actions/workflows/validate.yml)

Curated, provenance-tracked reference data for HLA–KIR–ERAP biological classifications used in immunogenetics and pharmacogenomics research.

## What this provides

Standardized lookup tables for NK cell education and antigen presentation pathway components, suitable for studies of drug hypersensitivity (DRESS, SJS/TEN), ankylosing spondylitis, viral immunity, and other HLA-mediated phenotypes.

### Database-derived classifications (managed)

| Dataset | Alleles | Source | Resolution |
|---------|---------|--------|------------|
| HLA-B Bw4/Bw6 (80I/80T) | 6,350 | IPD-IMGT/HLA 3.63.0 + HLAtools | 2-field |
| HLA-C C1/C2 | 4,978 | IPD-IMGT/HLA 3.63.0 + HLAtools | 2-field |
| HLA-B leader (-21M/T) | 4,353 | HLAtools protein alignment | 2-field |
| HLA-A + HLA-B Bw4/80I/80T | 7,953 | HLAtools protein alignment | 2-field |

### Literature-curated classifications

| Dataset | Entries | Primary source |
|---------|---------|----------------|
| KIR3DL1 expression (high/low/null) | 17 | Thomas 2008, Boudreau 2016 |
| KIR3DL1–HLA binding (full MFI matrix) | 323 | Maiers 2024 |
| ERAP1 allotype activity (10 allotypes) | 10 | Hutchinson 2021 |
| ERAP2 haplotype expression | 3 | Andrés 2010 |
| HLA functional divergence (pairwise FD) | 3,002 | Carrington 2024 |
| HLA-A estimated expression (z-scores) | 21 | Ramsuran 2018 |
| HLA-C allotype expression (MFI) | 14 | Petersdorf 2014 |
| ERAP SNP crosswalk (GRCh38 ↔ amino acid) | 11 | Ensembl REST VEP (r116) / dbSNP b156 |

## Quick start

### Use as a git submodule in your project

```bash
git submodule add https://github.com/krantzlab/immunogenomics-ref.git reference_data
```

### Load in R

```r
# Standalone:
source("R/load_reference_data.R")
# Or as submodule:
# source("reference_data/R/load_reference_data.R")

# Tidy tibbles for flexible analysis
bw4_data  <- load_bw4_classification()
c1c2_data <- load_c1_c2_classification()
kir_expr  <- load_kir3dl1_expression()
erap1     <- load_erap1_allotypes()
erap2     <- load_erap2_haplotypes()
hla_fd    <- load_hla_divergence()
hla_a_exp <- load_hla_a_expression()
hla_c_exp <- load_hla_c_expression()

# Or use backward-compatible vector extractors
BW4_80I_ALLELES <- get_bw4_80i_alleles()    # character vector of B* alleles
C1_ALLELES      <- get_c1_alleles()          # character vector of C* alleles
KIR3DL1_HIGH    <- get_kir3dl1_high()        # character vector of allele codes
B_LEADER        <- get_b_leader_lookup()     # named character vector (M/T)
```

### Validate data integrity

```r
source("R/validate_reference_data.R")
validate_all()
#   PASS: bw4_bw6
#   PASS: bw4_80i
#   PASS: c1_c2
#   PASS: b_leader
#   PASS: kir3dl1_expression
#   PASS: kir3dl1_binding
#   PASS: erap1_activity
#   PASS: erap2_expression
#   PASS: erap_crosswalk
#   PASS: hla_a_expression
#   PASS: hla_c_expression
#   PASS: hla_divergence
#   PASS: backward_compat
#   All reference data checks passed.
```

## Data architecture

```
immunogenomics-ref/
├── datasets/          # ← Downstream projects read ONLY from here
├── managed/           # R scripts: IPD-IMGT/HLA API → datasets/
├── curated/           # Literature CSVs with provenance → copy to datasets/
├── R/                 # Loader + validation functions
├── provenance.yaml    # Full source registry with DOIs and key decisions
└── CHANGELOG.md
```

**Managed data** is derived from the IPD-IMGT/HLA Allele Query API and cross-validated against HLAtools protein alignments. Update by running the fetch scripts when IPD releases quarterly updates.

**Curated data** is manually maintained from published literature. Every entry includes `confidence_tier` (A–D), `source_doi`, `source_detail`, and `evidence_note` for full traceability.

See `provenance.yaml` for the complete registry of sources, analytical decisions, and their rationale.

## Updating

### Refresh database-derived data

```bash
Rscript managed/fetch_kir_ligand.R    # Bw4/Bw6 + C1/C2
Rscript managed/fetch_b_leader.R      # -21 leader
Rscript managed/derive_bw4_80i.R      # Combined Bw4-80I/80T (A+B)
```

### Update curated data

1. Edit the CSV in `curated/`
2. Copy to `datasets/`
3. Update `provenance.yaml`
4. Run `validate_all()`
5. Update `CHANGELOG.md`

## R dependencies

**Core** (for loading and validation): `readr`, `dplyr`, `stringr`, `here`

**Managed scripts** (for fetching): additionally `httr`, `jsonlite`, `HLAtools`

### Development environment (pixi)

The consumer path needs none of this — `datasets/` is plain CSV, readable by
any language. `pixi.toml` exists so the loaders, validators and fetch scripts
can actually be *run*, with the interpreter and package versions pinned by
`pixi.lock`.

```bash
pixi install          # default environment: R 4.5 + readr, dplyr, stringr, here
pixi run validate     # all 13 checks; exits non-zero on failure
pixi run console      # interactive R with the loaders available
```

The `managed` environment adds `httr`/`jsonlite` for the IPD fetch scripts.
`HLAtools` is on CRAN but has no conda-forge build, so it installs into a
project-local `.Rlib/`:

```bash
pixi run -e managed deps      # one-time: install HLAtools into .Rlib/
pixi run -e managed refresh   # all three fetch scripts, then validate
```

`refresh` **rewrites files in `datasets/` in place** — that is the canonical
data every downstream consumer reads. Run it when IPD publishes a quarterly
release, and diff the result before committing.

## Key references

- Thomas R et al. (2008) *J Immunol* 180:6743. [KIR3DL1 expression](https://doi.org/10.4049/jimmunol.180.10.6743)
- Boudreau JE et al. (2016) *J Immunol* 196:189. [KIR3DL1 functional reclassification](https://doi.org/10.4049/jimmunol.1502469)
- Maiers M et al. (2024) *bioRxiv* / *J Biol Chem* (2025). [KIR3DL1–HLA binding quantification](https://doi.org/10.1101/2024.05.03.592082)
- Hutchinson JP et al. (2021) *J Biol Chem* 296:100443. [ERAP1 allotype activity](https://doi.org/10.1016/j.jbc.2021.100443)
- Andrés AM et al. (2010) *PLoS Genet* 6:e1001029. [ERAP2 balancing selection](https://doi.org/10.1371/journal.pgen.1001029)
- Foley BA et al. (2008) *J Immunol* 180:3969. [HLA-A Bw4 as KIR3DL1 ligand](https://doi.org/10.4049/jimmunol.180.6.3969)
- Petersdorf EW et al. (2022) *Blood* 139:3022. [HLA-B leader and NK education](https://doi.org/10.1182/blood.2021014437)
- Ramsuran V et al. (2018) *Science* 359:86. [HLA-A expression and HIV control](https://doi.org/10.1126/science.aam8825)
- Petersdorf EW et al. (2014) *Blood* 124:3996. [HLA-C expression levels in HCT](https://doi.org/10.1182/blood-2014-09-599969)
- Apps R et al. (2013) *Science* 340:87. [HLA-C expression and HIV control](https://doi.org/10.1126/science.1232685)
- Carrington M et al. (2024) *Science*. [HLA functional divergence](https://doi.org/10.1126/science.ado8609)

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on adding or updating reference data.

## License

Code and the compilation of these tables: [MIT](LICENSE).

The data itself carries the terms of its sources. Managed tables are derived
from IPD-IMGT/HLA; curated tables are extracted from published literature,
several from supplementary material under publisher copyright. Every curated
row carries `source_doi` so any individual value is traceable.

**See [DATA_TERMS.md](DATA_TERMS.md)** for the per-file breakdown. Short
version: cite these files freely, cite the original publication when a specific
value matters to your result, and check the source's terms before
re-publishing a curated table as your own dataset.

## Acknowledgments

- **IPD-IMGT/HLA Database** and the **EBI Allele Query API** for authoritative HLA classification data
- **HLAtools** R package (ANHIG/IMGTHLA) for protein sequence alignments
- **Claude** (Anthropic) — AI assistant used during development of the data architecture, validation framework, provenance system, and documentation for this repository
