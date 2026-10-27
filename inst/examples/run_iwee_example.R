#!/usr/bin/env Rscript
# Example: run IWEE() on a simulation file named like
#   paras_nTaxa<nTaxa>_<seed>_<Xtype>_X_<SNR>_<metric>[_<k>].rds
# Usage:   Rscript run_iwee_example.R [path/to/paras_file.rds] [num_cores]
# Without arguments the example data shipped with the package are used.
# This is the package version of run_iwee_diffzero_rep.R: it attaches IWEE
# instead of source()-ing the utility script and writes results to ./results.

library(IWEE)

args <- commandArgs(trailingOnly = TRUE)
rdsfile <- if (length(args) >= 1) args[1] else
  system.file("extdata", "iwee_example.rds", package = "IWEE")
num_cores <- if (length(args) >= 2) as.integer(args[2]) else 2L

## ---- tuning parameters by number of taxa (as in the original script) ----
nTaxa <- 400
randseed <- 1
if (grepl("paras_nTaxa", basename(rdsfile))) {
  parts <- regmatches(basename(rdsfile),
                      regexec("paras_nTaxa(\\d+)_(\\d+)_", basename(rdsfile)))[[1]]
  nTaxa <- as.numeric(parts[2])
  randseed <- as.numeric(parts[3])
}

settings <- list(
  `400` = list(nGroup = 16, Ng = 50, nRef = 10, nBestRef = 4),
  `200` = list(nGroup = 4,  Ng = 50, nRef = 40, nBestRef = 2),
  `50`  = list(nGroup = 2,  Ng = 24, nRef = 20, nBestRef = 2),
  `3`   = list(nGroup = 1,  Ng = 3,  nRef = 2,  nBestRef = 1)
)
cfg <- settings[[as.character(nTaxa)]]
if (is.null(cfg)) stop("No tuning parameters defined for nTaxa = ", nTaxa)

paras_obj <- readRDS(rdsfile)

start_time <- Sys.time()
iwee_res <- IWEE(
  paras_obj$MicrobData,
  paras_obj$CovData,
  num_cores = num_cores,
  testCov = "X1",
  ctrlCov = NULL,
  linkIDname = "id",
  QRdecomp = TRUE,
  nEta = 11,
  Ng = cfg$Ng,
  nGroup = cfg$nGroup,
  nRef = cfg$nRef,
  nBestRef = cfg$nBestRef,
  k_fold = 2,
  testMany = FALSE,
  ctrlMany = TRUE,
  randseed = randseed
)
runtime <- Sys.time() - start_time

## ---- summarize ----
cat("Runtime (secs):", as.numeric(runtime, units = "secs"), "\n")
cat("Taxa flagged as associated (FDR < 0.05):",
    sum(iwee_res$final_sig_beta), "\n")

if (!is.null(paras_obj$truBeta)) {
  truth <- paras_obj$truBeta["X1", names(iwee_res$final_beta_est)]
  cat("Taxa with a true X1 effect:", sum(truth != 0), "\n")
  cat("Of which detected:",
      sum(iwee_res$final_sig_beta & truth != 0), "\n")
}

dir.create("results", showWarnings = FALSE)
saveRDS(iwee_res, file.path("results", sub("^paras", "iwee", basename(rdsfile))))
