#!/usr/bin/env Rscript
#
# Run every reference-data check and exit non-zero if any of them fails.
#
# validate_all() prints its own PASS/FAIL lines, but it returns its result
# invisibly and never signals a condition. Calling it directly from a shell
# therefore always exits 0, and a failing dataset would show a green CI check.
# This wrapper is what makes the failure visible.
#
# Both CI (.github/workflows/validate.yml) and `pixi run validate` go through
# this file, so the two cannot drift apart.
#
# Usage:
#   Rscript tools/validate.R      # from the repository root

entry <- "R/validate_reference_data.R"

if (!file.exists(entry)) {
  stop(sprintf("%s not found. Run this from the repository root.", entry),
       call. = FALSE)
}

source(entry)

ok <- validate_all()

if (!isTRUE(ok)) {
  message("\nReference data validation FAILED -- see the FAIL lines above.")
  quit(status = 1L)
}

quit(status = 0L)
