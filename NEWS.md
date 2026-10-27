# IWEE 0.1.0

* First packaged release, converted from the script
  `IWEE_utils_withoutSparsematrix_QRdecomp_v4.R`.
* The functions are grouped in `R/IWEE.R` (main function), `R/parallel_jobs.R`
  (worker functions), `R/preprocessing.R` (`metaData()`, `estiEta()`) and
  `R/alpha_estimation.R` (estimation routines).
* `IWEE()`: `final_beta_pval`, `final_beta_fdr` and `final_sig_beta` now carry
  the original taxon names, like `final_beta_est` and `final_beta_var`.
* `IWEE()`: `num_cores` now defaults to `NA` (physical cores), as the code
  already assumed.
* Fixed a misspelled variable (`paris_count`) in the internal, unused
  `parallel_job_cv_v2()`.
* Added `load_iwee_example()` and a small simulated example dataset.
* Removed the `library()` calls and the unused `peakRAM` dependency; the
  commented-out Newton-Raphson code was not carried over.
