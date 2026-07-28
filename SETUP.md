# Environment setup

Tested on macOS 26.5.2 (arm64), R 4.6.1, via Homebrew R (`/opt/homebrew/bin/R`).

## Prerequisites

1. **Xcode license must be accepted** before R can compile any package from
   source. If `install.packages()` fails with errors like `had non-zero exit
   status` and the log mentions "You have not agreed to the Xcode license
   agreements", run:

   ```bash
   sudo xcodebuild -license
   ```

   and accept it (requires an admin password — this step can't be scripted
   around).

2. **cmake is required** to build `RcppParallel` (a dependency of
   `rxode2ll`/`rxode2`). Without it, install fails with `cmake was not found`.

   ```bash
   brew install cmake
   ```

## Installing R packages

```bash
Rscript -e 'options(repos=c(CRAN="https://cloud.r-project.org")); install.packages(c("rxode2","ggplot2"))'
```

**Do not pass `dependencies = TRUE`.** It recursively pulls in every
`Suggests` of every dependency — in practice this dragged in unrelated,
heavyweight packages (e.g. a full DuckDB source build) that have nothing to
do with running rxode2 models, and can take 10-20x longer to install. The
default `dependencies = NA` (Depends/Imports/LinkingTo only) is sufficient.

## Verifying the install

```bash
Rscript -e 'inst <- installed.packages()[,"Package"]; for (p in c("rxode2","ggplot2")) cat(p, ":", p %in% inst, "\n")'
```

## Known harmless warning

rxode2 prints this on every run on macOS; it does not affect correctness,
only single- vs multi-threaded solving speed:

```
rxode2 has not detected OpenMP support and will run in single-threaded mode
This is a Mac. Please read https://mac.r-project.org/openmp/
```
