# Data terms

This repository contains three kinds of material with different provenance, and
the [MIT licence](LICENSE) does not apply uniformly to all of them. This file
records what is licensed by whom.

It is a statement of provenance, not legal advice. If you intend to
redistribute any of the curated tables, check the terms of the underlying
publication.

## 1. Code — MIT

`R/`, `managed/`, `tools/`, `.github/`, and the build configuration are
original work and are licensed MIT, as stated in `LICENSE`.

## 2. Managed data — derived from IPD-IMGT/HLA

`datasets/bw4_bw6_classification.csv`, `datasets/c1_c2_classification.csv`,
`datasets/b_leader_assignments.csv`, `datasets/bw4_80i_classification.csv`

These are computed from the IPD-IMGT/HLA database, via its Allele Query API and
via protein alignments accessed through the `HLAtools` package. They are
derived facts — a classification computed per allele — rather than copies of
IPD records, but they exist only because IPD publishes the underlying data.

Users of these files should cite IPD-IMGT/HLA and observe its terms:

- Robinson J, Barker DJ, Georgiou X, Cooper MA, Flicek P, Marsh SGE.
  *IPD-IMGT/HLA Database.* Nucleic Acids Research (2020) 48:D948–D955.
  <https://doi.org/10.1093/nar/gkz950>
- <https://www.ebi.ac.uk/ipd/imgt/hla/licence/>

The `ipd_version` column in each file records the release the classifications
were derived from.

## 3. Curated data — extracted from published literature

Each row of every curated table carries `source_doi`, `source_detail` and
`confidence_tier`, so the origin of any individual value is traceable from the
data itself.

These values were extracted from tables in the following publications, several
of them from supplementary material. **Copyright in that supplementary material
is generally held by the publisher, not by this repository**, and the MIT
licence is not a grant over it.

| File | Rows | Source |
|---|---|---|
| `hla_functional_divergence.csv` | 3,002 | Carrington et al. (2024) *Science* — [10.1126/science.ado8609](https://doi.org/10.1126/science.ado8609) |
| `kir3dl1_hla_binding.csv` | 323 | Maiers et al. (2024) supplementary table — [10.1101/2024.05.03.592082](https://doi.org/10.1101/2024.05.03.592082) |
| `hla_a_estimated_expression.csv` | 21 | Ramsuran et al. (2018) *Science* Table S1 — [10.1126/science.aam8825](https://doi.org/10.1126/science.aam8825) |
| `hla_c_expression.csv` | 14 | Petersdorf et al. (2014) *Blood* Table 3 — [10.1182/blood-2014-09-599969](https://doi.org/10.1182/blood-2014-09-599969) |
| `kir3dl1_expression.csv` | 17 | Thomas et al. (2008) *J Immunol* Table II — [10.4049/jimmunol.180.10.6743](https://doi.org/10.4049/jimmunol.180.10.6743) |
| `erap1_allotype_activity.csv` | 10 | Hutchinson et al. (2021) *JBC* — [10.1016/j.jbc.2021.100443](https://doi.org/10.1016/j.jbc.2021.100443) |
| `erap2_haplotype_expression.csv` | 3 | Andrés et al. (2010) *PLoS Genet* — [10.1371/journal.pgen.1001029](https://doi.org/10.1371/journal.pgen.1001029) |
| `erap_snp_crosswalk.csv` | 11 | Ensembl REST VEP (release 116) over dbSNP build 156 |

**What this repository does claim:** the selection, arrangement, harmonisation,
provenance annotation and validation of these tables — the compilation — is
original work and is offered under MIT.

**What it does not claim:** the underlying measured values. Cite the original
publication, not this repository, when a specific value matters to your result.
`CITATION.cff` lists all of them under `references`.

## If you are redistributing

- Citing or analysing these files: no issue; cite the sources.
- Re-publishing a curated table as your own dataset: check the source
  publication's terms first. Several are behind publisher copyright.
- Vendoring a snapshot into a software package (as `bridgie` does): keep this
  file, or an equivalent statement, alongside the copied data.

## Corrections

If you hold rights in any material here and consider its inclusion improper,
open an issue and it will be removed promptly.
