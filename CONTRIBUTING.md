# Contributing to immunogenomics-ref

Thank you for your interest in contributing. This repository maintains reference data that must be accurate and well-sourced, so contributions follow a structured process.

## Types of contributions

### 1. Adding new curated datasets

New literature-sourced classification tables are welcome. Each must include:

- **CSV file** in `curated/` with required provenance columns: `confidence_tier`, `source_doi`, `source_detail`, `evidence_note`
- **Loader function** (`load_*()`) in `R/load_reference_data.R` with column validation
- **Validation function** (`validate_*()`) in `R/validate_reference_data.R` with spot-checks
- **Provenance entry** in `provenance.yaml` with full source details, analytical decisions, and confidence ratings
- **Getter functions** (`get_*()`) if backward-compatible vector extractors are needed
- **CHANGELOG.md** entry

### 2. Updating existing classifications

When new literature changes a classification:

1. Edit the CSV in `curated/`
2. Update `provenance.yaml` — add the new reference, update `last_literature_review`, document the decision change
3. Copy updated CSV to `datasets/`
4. Run `validate_all()` — all checks must pass
5. Update `CHANGELOG.md`

Please include the DOI and a brief rationale in your pull request description.

### 3. Refreshing managed data

When IPD-IMGT/HLA releases a new version:

1. Run the fetch scripts in order
2. Note any new cross-validation mismatches
3. Update the "Last fetch" section in `managed/README.md`
4. Run `validate_all()`
5. Update `CHANGELOG.md`

### 4. Bug fixes and improvements

Issues with validation logic, loader functions, or documentation are welcome as PRs.

## Confidence tier guidelines

When classifying evidence for curated entries:

| Tier | Use when... |
|------|-------------|
| A | You have direct experimental measurement for this specific allele/allotype from the cited study |
| B | Multiple studies agree, or the classification is well-established in the field |
| C | Classification is inferred from a related allele's data (e.g., subfamily similarity), or the data is limited |
| D | Single unreplicated finding, or expert opinion without direct experimental support |

## Commit conventions

Use conventional commits:
- `data:` — changes to dataset files
- `data(curated):` — curated data updates (include DOI)
- `data(managed):` — managed data refreshes (include IPD version)
- `feat:` — new loader/validator functions
- `fix:` — bug fixes
- `docs:` — documentation changes
- `build:` — development environment and tooling (`pixi.toml`, `tools/`)
- `ci:` — GitHub Actions workflows

## Pull request checklist

- [ ] `pixi run validate` passes (13/13 checks)
- [ ] New curated data includes all provenance columns
- [ ] `provenance.yaml` updated for any data changes
- [ ] `CHANGELOG.md` updated
- [ ] DOIs are valid and resolve correctly
- [ ] No study-specific data included (this repo is study-agnostic)
- [ ] No AI attribution trailers (`Co-authored-by`, `Claude-Session:`, "Generated with")

## Attribution

AI assistants are not attributed as commit co-authors. Git attribution implies
accountability for the committed content, and an assistant cannot hold it;
assistant contributions are acknowledged in `README.md` instead. Enable the
enforcing hook once per clone with `git config core.hooksPath .githooks`. If it
rejects a commit, remove the trailer rather than bypassing with `--no-verify`.

## Release checklist

Downstream packages pin this repository by tag and vendor snapshots of
`datasets/` (for example `bridgie` reads a pinned copy from its own
`inst/extdata/`). A published tag is therefore part of their provenance record
and must never be moved, retagged or deleted — only superseded by a new tag.

- [ ] `pixi run validate` passes
- [ ] `CHANGELOG.md` has a `## vX.Y.Z` section for the release
- [ ] `CITATION.cff` `version` and `date-released` match the tag
- [ ] Schema changes are reflected in the version bump: additive columns are a
      minor release; changing a column's value domain or renaming/removing a
      column is breaking and needs a major release plus a `CHANGELOG.md` note
      naming the affected consumers
- [ ] Tag is annotated (`git tag -a vX.Y.Z -m ...`) and pushed with
      `git push origin vX.Y.Z`
- [ ] The release notes quote the **commit** SHA, not the annotated tag object
      SHA, so consumers record an unambiguous pin
