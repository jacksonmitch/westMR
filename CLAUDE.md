# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Current focus

Only working on documentation right now: roxygen2 comments, README, 
vignettes, NEWS.md. Do not modify function logic in R/ unless explicitly asked.

## Git workflow

- Never commit directly to main
- Always work on a feature branch, open a PR when done
- Don't use git add -A; stage files selectively
- This repo is shared with a research partner — keep changes scoped and 
  easy to review

## What this package does

westMR is an R package implementing the WEST procedure for finite mixture
regression (FMR): stepwise selection of which predictors belong in the model,
and of whether each predictor's effect is homogeneous (shared across mixture
components) or heterogeneous (component-specific). Supports Gaussian,
Poisson, and binomial response families. Estimation uses an EM algorithm,
with weighted least squares (via structured QR decomposition) in the M-step,
and IRWLS for non-Gaussian families.

## Common commands

Run these from an R console with the working directory set to the package
root (`devtools`/`testthat` are expected to be available via `renv`).

```r
# Load package for interactive development
devtools::load_all()

# Run the full test suite
devtools::test()

# Run a single test file
testthat::test_file("tests/testthat/test-determine_effects.R")

# Run tests matching a name/description
devtools::test(filter = "determine_effects")

# Regenerate NAMESPACE and man/ pages from roxygen comments (after editing
# any roxygen block, e.g. in westMR.R or control.R)
devtools::document()

# Full package check (mirrors the R-CMD-check.yaml CI workflow)
devtools::check()

# Lint (mirrors the lint.yaml CI workflow; must pass with zero lints)
lintr::lint_package()
```

renv manages dependencies (`renv.lock`); `.Rprofile` activates it
automatically when R starts in this project.

CI (`.github/workflows/`) runs `R-CMD-check.yaml` (check across
macOS/Windows/Ubuntu, release/devel/oldrel) and `lint.yaml`
(`lintr::lint_package()`, errors on any lint) on every push/PR to
main/master.

## Architecture

### Entry point and flow

`westMR()` (R/westMR.R) is the public entry point. It validates arguments,
optionally coerces a binomial formula (`coerce_binomial_formula`, handles
both `cbind(success, failure)` and 0/1 responses by attaching a
`.binom_size` column to the data), builds a `WMRModel` (R6, R/WMRModel.R),
then runs up to two stepwise procedures depending on `task`:

1. `select_variables()` (R/select_variables.R) — decides which predictors
   are in vs. out of the model.
2. `determine_effects()` (R/determine_effects.R) — decides, for included
   predictors, which are heterogeneous vs. homogeneous.

Both are thin wrappers around the shared driver `test_predictors()`
(R/test_predictors.R), which performs greedy forward/backward stepwise
search: at each step it fits a "shared" model plus one candidate model per
remaining predictor, runs `west_procedure()` to compare each candidate
against the shared fit via a BIC-weighted likelihood ratio test, and either
adopts the most eligible candidate or stops. `select_variables` and
`determine_effects` differ only in what changes at each step (included set
vs. heterogeneous set) — this is parameterized via the `update_included` /
`heterogeneous` arguments passed into `test_predictors`.

`west_procedure()` (R/west_procedure.R) computes, for every G in
`model$G_values`, a null-vs-alternative BIC/LRT comparison, then combines
the per-G p-values into a single BIC-weighted p-value (`p0`) used for the
stepwise eligibility decision (`p0 < alpha` forward / `p0 >= alpha`
backward).

### Model fitting stack (bottom-up)

- **`WMRData`** (R/prepare_data.R, R6) — validated container for one
  design-matrix specification: `X_het` (heterogeneous/component-specific
  design matrix) and `X_com` (common/homogeneous design matrix), built by
  `prepare_data(model, included, common)` from the model formula. Read-only
  fields prevent accidental mutation after construction.
- **`EmState`** (R/em_state.R, R6) — mutable container for one EM run's
  state (`tau`, `pi_g`, `beta_g`, `beta`, `sigma_g`, `eta`, `loglik`, ...).
  Passed by reference into `m_step`/`e_step_fmr`/`irwls_fmr`, which mutate
  it in place rather than returning a new object each iteration.
- **`em_fmr()`** (R/EM_fmr.R) — the EM loop: alternates `m_step()` and
  `e_step_fmr()` until the log-likelihood change falls below `control$tol`
  or `control$max_iter` is reached.
  - **`m_step()`** (R/M_step.R) dispatches on family, and on whether any
    common/homogeneous predictors exist (`dat$p_com > 0`): with common
    predictors it uses the structured-QR path (`m_step_sqr_*`,
    `wls_sqr()` in R/wls_sqr.R, which jointly solves for heterogeneous and
    common coefficients component-by-component using QR); without common
    predictors it falls back to plain per-component `lm.wfit`/`glm.fit`.
    Poisson/binomial M-steps run an inner IRWLS loop (`irwls_fmr()` in
    R/irwls.R) that repeatedly re-weights and re-solves via `wls_sqr()`.
  - **`e_step_fmr()`** (R/E_step.R) computes per-component log-densities
    (family-specific: dnorm/dpois/dbinom), updates responsibilities `tau`
    via log-sum-exp (numerically stable, R/utils.R `row_logsumexp`), and
    updates mixing proportions `pi_g`.
- **`fit_fmr()`** (R/fit_fmr.R) fits one model spec at one G: it first
  burns in several `tau` initializations (`make_tau_list()` in
  R/initialization.R — quantile-based on response/residuals, k-means, plus
  random balanced partitions; for binomial also per-predictor quantile
  splits) for a few EM iterations each via `select_best_initialization()`,
  picks the best by log-likelihood, then runs `em_fmr()` to convergence
  from that state. Also computes `k` (parameter count, `count_params_gmr()`
  in R/model_criteria.R) and BIC (`compute_bic()`).
  `fit_across_G()` calls `fit_fmr()` once per value in `model$G_values` and
  optionally seeds one extra initialization per G from a previously-fit
  `tau` (`extra_tau_starts`, used to warm-start candidate fits from the
  shared fit in `west_procedure`).

### Class/return-object conventions

Every major result carries an S3 class (`fit_fmr`, `select_variables`,
`determine_effects`, `westMR`, `west_procedure`) with print/summary methods
defined centrally in R/s3-methods.R, not alongside the constructors. The
top-level `westMR` object's `print`/`summary` recursively delegate to its
components' own print methods. `steps_table()` is a generic (declared in
R/s3-methods.R) for rendering the step-by-step stepwise search log stored on
`select_variables`/`determine_effects` objects.

R6 classes (`WMRModel`, `WMRData`, `EmState`) are used where mutable
reference semantics or derived/read-only active bindings are needed (e.g.
`WMRModel$predictors`/`response` are derived from `formula` and cannot be
set directly). Plain S3 lists are used for immutable result objects.

### Numerics notes worth knowing before editing

- `wls_sqr()` implements weighted least squares via a two-stage QR
  (component-specific QR to eliminate heterogeneous coefficients, then a
  pooled QR to solve for common coefficients) rather than forming normal
  equations directly — this is the core "structured QR" referenced in
  DESCRIPTION.
- Numerical floors (`sigma_floor`, `weight_floor`, `1e-16`/`1e-8` clamps)
  appear throughout the E/M steps to avoid degenerate components; when
  touching these paths, preserve the existing floor pattern rather than
  removing clamps.
- `control$sigma_floor`, when left `NULL` in `build_control()`, is computed
  lazily in `WMRModel$initialize()` as `0.05 * sd(response)`.

### Binomial family

The binomial path expects a `.binom_size` column added to `data` by
`coerce_binomial_formula()` (called from `westMR()` before `WMRModel` is
built); `.binom_size` is a reserved column name. `m_step_binomial()` (plain,
non-QR path) is not implemented (`stop("binomial not done yet")`) — only
the structured-QR + IRWLS path (`m_step_sqr_binomial`) currently works.

## Tests

`tests/testthat/helper.R` defines `simulate_fmr()` and a `scenarios` list
(pre-defined mixture configurations: two/three/four-group, subtle vs.
strong separation, unbalanced mixing proportions, heterogeneous variance)
used across test files to generate synthetic FMR data with known ground
truth. Reuse these scenarios for new tests rather than inventing new
simulation setups. Testthat edition 3 (`Config/testthat/edition: 3` in
DESCRIPTION).
