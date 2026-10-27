test_that("IWEE() runs end to end and returns named results", {
  skip_on_cran()
  dat <- load_iwee_example()
  n_taxa <- 60
  out <- capture.output(
    res <- IWEE(
      MicrobData = dat$MicrobData[, seq_len(n_taxa + 1)],
      CovData = dat$CovData,
      num_cores = 2,
      testCov = "X1",
      ctrlCov = c("W1", "W2"),
      linkIDname = "id",
      QRdecomp = TRUE,
      nEta = 3,
      Ng = 30,
      nRef = 5,
      nBestRef = 2,
      k_fold = 2,
      testMany = FALSE,
      ctrlMany = TRUE,
      randseed = 1
    )
  )
  expected <- c("eta_cvMSE", "etaMSE", "etaCand", "bestEtaCount", "bestEtaInd",
                "final_ref_names.lst", "final_beta_var", "final_beta_est",
                "final_beta_pval", "final_beta_fdr", "final_sig_beta",
                "randseed", "phase1_count_coefMagnit")
  expect_true(all(expected %in% names(res)))
  expect_length(res$etaCand, 3)
  expect_equal(res$randseed, 1)

  # all result vectors are named with the original taxon names
  est_names <- names(res$final_beta_est)
  expect_true(all(est_names %in% paste0("taxon", 1:400)))
  expect_equal(names(res$final_beta_pval), est_names)
  expect_equal(names(res$final_beta_fdr), est_names)
  expect_equal(names(res$final_sig_beta), est_names)
  expect_type(res$final_sig_beta, "logical")
})
