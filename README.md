# IWEE

Inverse-weighted estimating equations (IWEE) for identifying taxa associated
with covariates in zero-inflated, compositional microbiome data.

## Installation

From the folder that contains the `IWEE` directory:

```r
install.packages(c("doParallel", "foreach", "glmnet", "MASS", "Matrix"))
install.packages("IWEE", repos = NULL, type = "source")
# or: R CMD INSTALL IWEE
```

Alternatively, install directly from GitHub:

```r
install.packages("remotes")
remotes::install_github("TChang-W/IWEE-microbiome")
```

The parallel stages start PSOCK workers that load the *installed* package, so
install it (do not only `devtools::load_all()`) before calling `IWEE()`.

## Example

```r
library(IWEE)

dat <- load_iwee_example()   # 100 samples, 400 taxa, covariates X1, W1, W2

# Reduced run (first 100 taxa, small tuning grid)
res <- IWEE(
  MicrobData = dat$MicrobData[, 1:101],
  CovData    = dat$CovData,
  num_cores  = 2,
  testCov    = "X1",
  ctrlCov    = c("W1", "W2"),
  linkIDname = "id",
  QRdecomp   = TRUE,
  nEta       = 3,
  Ng         = 50,
  nRef       = 10,
  nBestRef   = 2,
  k_fold     = 2,
  randseed   = 1
)

head(res$final_beta_est)
names(which(res$final_sig_beta))   # taxa with FDR < 0.05
```

The full-size analysis of the example data (16 groups of 50 taxa) is in
`system.file("examples", "run_iwee_example.R", package = "IWEE")`.

## Contents

| File | Content |
|------|---------|
| `R/IWEE.R` | `IWEE()`, the main function |
| `R/parallel_jobs.R` | worker functions run by the parallel stages |
| `R/preprocessing.R` | `metaData()` (data cleaning), `estiEta()` |
| `R/alpha_estimation.R` | OLS/QR and LASSO estimation routines |
| `R/data.R` | `load_iwee_example()` |
