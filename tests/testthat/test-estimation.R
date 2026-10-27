make_toy_fit_data <- function(n = 200, seed = 3) {
  set.seed(seed)
  micro <- matrix(rexp(n * 4, rate = 0.5) + 0.1, n, 4,
                  dimnames = list(NULL, paste0("microb", 1:4)))
  cov <- cbind(1, x1 = rnorm(n))
  list(micro = micro, cov = cov)
}

test_that("InfCorrect() replaces infinite values and keeps the sign", {
  expect_equal(InfCorrect(c(1, Inf, -Inf)), c(1, 1e30, -1e30))
  expect_equal(InfCorrect(c(2, 3)), c(2, 3))
})

test_that("alpha_QRdecompose() solves a well-conditioned linear system", {
  A <- matrix(c(2, 1, 1, 3), 2)
  y <- c(1, 2)
  expect_equal(as.numeric(alpha_QRdecompose(A, y)), as.numeric(solve(A, y)),
               tolerance = 1e-6)
})

test_that("alpha_phi_Update_iloop_OLS() returns finite estimates for both solvers", {
  d <- make_toy_fit_data()
  for (qr_flag in c(TRUE, FALSE)) {
    res <- alpha_phi_Update_iloop_OLS(
      gMicrobData = d$micro, CovDataWithIntcp = d$cov, eta = 0.5,
      RefTaxonName = "microb4", QRdecomp = qr_flag
    )
    expect_length(res$alpha.hat, ncol(d$cov) * 3)
    expect_true(all(is.finite(res$alpha.hat)))
    expect_equal(dim(res$phiMat), c(4, 4))
    expect_equal(res$twoPosSub, seq_len(nrow(d$micro)))
  }
})

test_that("alpha_phi_Update_iloop_LASSO() returns one coefficient per covariate and pair", {
  d <- make_toy_fit_data()
  res <- alpha_phi_Update_iloop_LASSO(
    gMicrobData = d$micro, CovDataWithIntcp = d$cov, eta = 0.5,
    RefTaxonName = "microb4", testCovInd = 2
  )
  expect_length(res$alpha.hat, ncol(d$cov) * 3)
})
