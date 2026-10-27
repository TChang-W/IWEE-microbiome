make_toy_data <- function(n = 30, p = 6, seed = 1) {
  set.seed(seed)
  id <- seq_len(n)
  micro <- matrix(rpois(n * p, 3), n, p,
                  dimnames = list(NULL, paste0("t", seq_len(p))))
  micro[sample(length(micro), 20)] <- 0
  micro[1, ] <- seq_len(p) # no taxon or subject without reads
  cov <- cbind(id = id, X1 = rep(0:1, length.out = n), W1 = rnorm(n))
  list(micro = cbind(id = id, micro), cov = cov)
}

test_that("metaData() cleans, renames and orders the data", {
  toy <- make_toy_data()
  out <- capture.output(
    meta <- metaData(toy$micro, toy$cov, linkIDname = "id",
                     testCov = "X1", ctrlCov = "W1",
                     testMany = FALSE, ctrlMany = FALSE)
  )
  expect_equal(meta$Mprefix, "microb")
  expect_equal(meta$covsPrefix, "x")
  expect_equal(meta$microbName, paste0("t", 1:6))
  expect_equal(meta$newMicrobNames, paste0("microb", 1:6))
  expect_equal(nrow(meta$data), 30)
  expect_true(all(c("id", paste0("microb", 1:6), "x1", "x2") %in% colnames(meta$data)))
  # the binary covariate is moved after the non-binary one
  expect_equal(meta$xNames, c("W1", "X1"))
  expect_equal(meta$testCovInNewNam, "x2")
})

test_that("metaData() validates its input", {
  toy <- make_toy_data()
  capture.output({
    expect_error(
      metaData(toy$micro, toy$cov, linkIDname = "id", testCov = "nope",
               testMany = FALSE, ctrlMany = FALSE),
      "not available"
    )
    expect_error(
      metaData(toy$micro, toy$cov, linkIDname = "id", testCov = NULL,
               testMany = FALSE, ctrlMany = FALSE),
      "No covariates"
    )
  })
})

test_that("estiEta() returns a positive root that decreases in the mechanism-zero share", {
  set.seed(2)
  x <- matrix(rpois(200 * 5, 1), 200, 5)
  eta_low <- estiEta(x, mechZeroPerc = 0.2)
  eta_high <- estiEta(x, mechZeroPerc = 0.8)
  expect_true(is.numeric(eta_low) && length(eta_low) == 1 && eta_low > 0)
  expect_true(eta_high > 0)
  expect_gt(eta_low, eta_high)
})

test_that("load_iwee_example() returns the documented structure", {
  dat <- load_iwee_example()
  expect_type(dat, "list")
  expect_equal(dim(dat$MicrobData), c(100, 401))
  expect_equal(colnames(dat$CovData), c("id", "X1", "W1", "W2"))
  expect_equal(dim(dat$truBeta), c(4, 400))
})
