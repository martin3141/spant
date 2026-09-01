# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

spant ("Spectroscopy Analysis Tools") is an R package for reading, visualising, and
processing Magnetic Resonance Spectroscopy (MRS) data. It implements fully automated
metabolite fitting (ABfit), retrospective frequency/phase correction (RATS), NIfTI-MRS
I/O, quantum-mechanical basis-set simulation, and voxel registration for partial-volume
correction. Standard R package layout: `R/` (source), `src/` (Fortran/C), `tests/testthat/`
(unit tests), `vignettes/` (Rmd tutorials), `man/` (generated Rd docs), `inst/extdata`
(example data files used in examples/tests).

## Common commands

Run these from an R console with the working directory set to the package root (or via
`Rscript -e '...'`). `devtools` and `roxygen2` are the standard tools used here.

- Load package for interactive development: `devtools::load_all()`
- Run the full test suite: `devtools::test()`
- Run a single test file: `devtools::test_active_file("tests/testthat/test_fitting.R")`
  or `testthat::test_file("tests/testthat/test_fitting.R")`
- Regenerate documentation (Rd files in `man/`, `NAMESPACE`) from roxygen2 comments:
  `devtools::document()` — always run this after adding/changing `#' @export` or any
  roxygen block, and commit the resulting changes under `man/` and `NAMESPACE`.
- Full CRAN-style check (required before a PR is considered mergeable):
  `devtools::check(args = c('--as-cran'))`
- Build the package: `devtools::build()`
- Compile Fortran/C sources in `src/` (only needed after editing `src/*.f`/`src/*.c`):
  `devtools::load_all()` triggers a recompile automatically, or run `R CMD SHLIB` manually.
- Render README from README.Rmd after editing it: `devtools::build_readme()` (README.md
  is generated — always edit `README.Rmd`, never `README.md` directly).

CI (`.github/workflows/check-standard.yaml`) runs `R-CMD-check` across macOS, Windows,
and Ubuntu (devel/release/oldrel) on pushes/PRs to `main`/`master`.

## Architecture

### The `mrs_data` object

The central data type is the S3 class `mrs_data`, constructed by the internal
`mrs_data()` function in `R/mrs_data_io.R`. It's a list with fields:
`data` (a 7-dimensional complex array ordered `dummy, X, Y, Z, dynamic, coil, FID`),
`ft` (transmitter frequency), `resolution` (7-element vector, last element is `1/fs`),
`ref` (ppm reference), `nuc` (resonant nucleus), `freq_domain` (per-dimension logical
flags for time vs frequency domain), `affine` (voxel-to-world transform, may be NULL),
`meta`, and `extra` (optional data frame for grouping/id variables in batch analyses).
This 7D array layout is what allows single-voxel, dynamic (time-series), and MRSI
(imaging) data to share the same representation and generic functions.

Helper constructors build `mrs_data` objects from simpler inputs without needing every
field: `vec2mrs_data`, `mat2mrs_data`, `array2mrs_data` (all in `R/mrs_data_proc.R`).
Validate an incoming object with `check_mrs_data()` rather than re-deriving the check.

Most processing/plotting functions are S3 generics with an `.mrs_data` method (e.g.
`lb.mrs_data`, `zf.mrs_data`, `ppm.mrs_data`, `plot.mrs_data`, arithmetic operators),
declared in `R/mrs_data_proc.R` (~5,600 lines — the largest file, containing most
pre-processing steps: phasing, coil combination, zero-filling, apodization, HSVD
filtering, alignment, etc). When adding a new processing step, follow this generic +
`.mrs_data` method pattern and register it via roxygen `@export`.

### I/O

Format-specific readers live in one file per vendor format:
`mrs_read_dicom.R`, `mrs_read_ima.R`, `mrs_read_jmrui_txt.R`, `mrs_read_lcm_raw.R`,
`mrs_read_list_data.R`, `mrs_read_nifti.R`, `mrs_read_paravis.R`, `mrs_read_pfile.R`
(GE), `mrs_read_rda.R` (Siemens), `mrs_read_spar_sdat.R` (Philips), `mrs_read_twix.R`
(Siemens raw), `mrs_read_varian.R`. All are dispatched through the single public
entry point `read_mrs()` in `R/mrs_data_io.R`, which guesses the format from the file
extension (or directory contents) unless `format` is given explicitly. Writers
(`mrs_write_nifti.R`) follow the inverse pattern. `dicom_reader.R` provides shared
low-level DICOM parsing used by several of the above.

### Fitting pipeline

- `R/fitting.R`: `fit_mrs()` is the top-level entry point, dispatching to a fitting
  `method` (default `"ABFIT"`), producing an S3 `fit_result` object (see
  `R/fit_display.R` for its structure and `plot.fit_result`/`stackplot.fit_result`).
- `R/abfit.R`: the ABfit algorithm — automated baseline + metabolite fitting.
- `R/varpro.R`, `R/varpro_basic.R`, `R/varpro_3_para.R`: VARPRO-family fitting methods.
- `R/tdsr.R`, `R/rats.R`: time-domain spectral registration / RATS frequency-and-phase
  alignment.
- `R/pnnls.R` wraps Fortran non-negative least squares solvers (`src/lsei.f`,
  `src/pnnls.f`, called via `.Fortran("nnls", ...)` / `.Fortran("pnnls", ...)`) used
  as fitting building blocks.
- `R/basis_set.R`: `basis_set` S3 class representing simulated/measured metabolite
  basis spectra used by the fitters; includes TARQUIN-based basis generation
  (`write_basis_tqn`, requires an external TARQUIN binary — see `check_tqn()` in
  `R/spant.R`) and quantum-mechanical simulation (`R/qm_simulation.R`,
  `R/pulse_sequences.R`, `R/pulse_shapes.R`, `R/gsl_functions.R`).
- `R/svs_batch_fit.R`: batch-processing pipeline for fitting many single-voxel datasets.
- `R/mol_parameters.R`: hard-coded chemical shift/coupling/relaxation parameters for
  standard brain metabolites, macromolecules, and lipids used in simulation and fitting.

### External tool integration

Some functionality shells out to external binaries configured via `options()`
(set in `.onLoad()` in `R/spant.R`): `spant.tqn_cmd` (TARQUIN, default `"tarquin"`)
and `spant.lcm_cmd` (LCModel). `check_tqn()`/`check_lcm()` verify these are runnable.
`R/image_reg.R` and `R/ants_utils.R` integrate with `RNiftyReg`/`rpyANTs` for voxel
registration to anatomical images (partial volume correction).

### Other notable files

- `R/precomp.R`: precompiled/cached basis or lookup data handling
  (`spant.precomp_mode`/`spant.precomp_verbose` options).
- `R/cli_utils.R`, `R/interactive_plotting.R`: CLI progress helpers and a Shiny-based
  interactive viewer (`Suggests`: shiny, shinyFiles, miniUI).
- `R/benchmark.R`: performance benchmarking helpers.
- `R/lofdc.R`: lipid/other signal removal or decomposition utility.
- `R/fmrs.R`: functional MRS (dynamic time-series) specific helpers.

## Testing

Tests use `testthat` (`tests/testthat.R` runs `test_check("spant")`). Test files are
grouped by area: `test_fitting.R`, `test_nifti.R`, `test_preproc.R`, `test_qm_sim.R`.
Several tests compare against cached `.rds` fixtures in `tests/testthat/` (e.g.
`abfit_res_default.rds`, `def_basis.rds`) — if you change fitting numerics, expect to
need to regenerate these fixtures deliberately, not incidentally. Example data used by
`@examples` blocks and tests lives under `inst/extdata`.

## Style

Follow the Hadley Wickham R style guide (`http://adv-r.had.co.nz/Style.html`), per
CONTRIBUTING.md. Roxygen2 with markdown support is used for all documentation
(`Roxygen: list(markdown = TRUE, roclets = c("rd", "namespace", "collate"))` in
`DESCRIPTION`) — never hand-edit `NAMESPACE` or files under `man/`; regenerate them
with `devtools::document()` instead.

## Versioning

`DESCRIPTION` `Version` and `Date` are bumped as part of the release process (see
recent commits like "version bump"); `NEWS.md` tracks user-facing changes per version.
