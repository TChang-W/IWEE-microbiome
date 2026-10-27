#' Replace infinite values by a large finite value
#'
#' Helper used to keep weighted cross-products finite: every infinite
#' element of `para` is replaced by `correct_scale` carrying its sign.
#'
#' @param para Numeric vector or matrix.
#' @param correct_scale Magnitude used in place of `Inf`/`-Inf`.
#'
#' @return `para` with infinite entries replaced.
#' @keywords internal
InfCorrect <- function(para, correct_scale = 1e+30) {
  para[which(is.infinite(para))] <- correct_scale*sign(para[which(is.infinite(para))])
  return(para)
}

#' Solve the normal equations with a QR decomposition
#'
#' Solves \eqn{A x = y} for \eqn{x}, where \eqn{A} is the cross-product matrix
#' `ZZsum` and \eqn{y} is `ZUsum`, using a sparse QR decomposition from
#' \pkg{Matrix}. If `ZZsum` is numerically rank deficient (an absolute
#' eigenvalue below `1e-10`), a ridge of `1e-10` is added to its diagonal
#' first.
#'
#' @param ZZsum Square cross-product matrix \eqn{\sum Z' W Z}.
#' @param ZUsum Cross-product vector or one-column matrix \eqn{\sum Z' W U}.
#' @param small_value Unused; retained from an earlier implementation.
#'
#' @return The coefficient estimate \eqn{x} as returned by `qr.coef()` on the
#'   sparse QR decomposition.
#' @keywords internal
alpha_QRdecompose <- function(ZZsum, ZUsum, small_value = 1e-10) {
  A <- Matrix(ZZsum, sparse = TRUE)
  if (any(abs(eigen(A, only.values = TRUE)$values) < 1e-10)) { # check if ZZsum is rank deficient
    A <- A + 1e-10 * diag(nrow(A))
  }
  qr.A <- Matrix::qr(A, tol = 0)
  # Q <- qr.Q(qr.A)
  # R <- qr.R(qr.A)
  # # if R is not invertible:
  # if (any(diag(R) == 0)) {
  #   return(NA)  # Return NA if R is singular
  # }
  # diag_R <- diag(R)
  # diag_R[diag_R == 0] <- small_value # I suspect that replacing 0 with non-zero value introduces noise to our methodology because of its poor performance of choosing best eta
  # diag(R) <- diag_R
  # 
  # Qy <- t(Q)%*%ZUsum
  # alpha.hat <- backsolve(R, Qy)
  alpha.hat <- qr.coef(qr.A, ZUsum)
  
  return(alpha.hat) # somtimes overfitting: max(alpha.hat)==69799621749: will ruin JK
}

#' Iterative weighted least squares for the IWEE coefficients
#'
#' Estimates the pairwise log-ratio regression coefficients \eqn{\alpha}
#' (covariate effects relative to the reference taxon) and the pair-specific
#' dispersions \eqn{\phi} by alternating between a weighted least-squares
#' update of \eqn{\alpha} and a moment update of \eqn{\phi}. Only subjects
#' with at least two positive taxa contribute. Weights are the inverse
#' probabilities of observing both taxa of a pair, determined by `eta`. The
#' loop stops after 10 iterations, once the mean absolute change in
#' \eqn{\alpha} is at most `1e-3`, or when the update is deemed abnormal.
#'
#' @param gMicrobData Numeric matrix (N by K+1) of taxa abundances with
#'   column names.
#' @param CovDataWithIntcp Numeric matrix (N by 1+Q) of covariates including
#'   an intercept column.
#' @param eta Parameter of the zero-generating mechanism (see [estiEta()]).
#' @param RefTaxonName Column name of the reference taxon in `gMicrobData`.
#' @param phiMat.init Initial (K+1) by (K+1) dispersion matrix.
#' @param printloss Logical. Reserved for diagnostics; must be `FALSE`.
#' @param abnormal_ginv_threshold When `QRdecomp = FALSE`, the update stops
#'   if the pseudo-inverse deviates from an inverse by more than this relative
#'   amount.
#' @param abnormal_QR_threshold When `QRdecomp = TRUE`, the update stops if
#'   the sum of absolute coefficients grows by more than this factor.
#' @param phi_control Pairs observed in at most `ncol(CovDataWithIntcp) +
#'   phi_control` subjects get no own dispersion estimate; they receive the
#'   mean of the other estimates.
#' @param QRdecomp Logical; solve with a QR decomposition (`TRUE`) or
#'   [MASS::ginv()] (`FALSE`).
#' @param jackknife Logical. If `TRUE`, the per-subject contributions to the
#'   cross-products are returned for use by the jackknife.
#'
#' @return A list with the per-subject lists `k2k1.lst`, `AZ.lst`, `W.lst`,
#'   `U.lst`, the subjects used (`twoPosSub`), the final `phiMat`, the
#'   cross-products `ZZsum` and `ZUsum`, the coefficient vector `alpha.hat`,
#'   the iteration histories `alpha_update_mat` and `phi_update_mat`, and,
#'   if `jackknife = TRUE`, the per-subject cross-products `AZWAZ.lst` and
#'   `AZWU.lst`. If no iteration is accepted, `alpha.hat` and `phiMat` are
#'   absent and only the initial estimate in `alpha_update_mat` is available.
#' @keywords internal
alpha_phi_Update_iloop_OLS <- function(
    gMicrobData, # N by (K+1) taxa matrix
    CovDataWithIntcp, # N by (1+Q+S) covariate matrix which includes intercept
    eta, # parameter of zero probability mechanism
    RefTaxonName, 
    phiMat.init = matrix(1, nrow = ncol(gMicrobData), ncol = ncol(gMicrobData)),
    printloss = FALSE,
    abnormal_ginv_threshold = 0.025,
    abnormal_QR_threshold = 10,
    phi_control = 2,
    QRdecomp = TRUE,
    jackknife = FALSE
){
  results <- list()
  force(phiMat.init)  # Ensures evaluation before parallel execution
  
  # calculate K: K+1 is the number of taxa
  K <- ncol(gMicrobData) - 1
  
  # reorganize the microbiome matrix by moving the reference taxon to the last column
  refTaxonInd <- which(colnames(gMicrobData)==RefTaxonName)
  org.gMicrobData <- cbind(gMicrobData[,-refTaxonInd], gMicrobData[,refTaxonInd])
  rm(gMicrobData)
  
  k2k1.lst <- list()
  AZ.lst <- list()
  W.lst <- list()
  U.lst <- list()
  ZZadd <- list()
  ZUadd <- list()
  twoPos.vec <- rep(FALSE, nrow(org.gMicrobData))
  for (i in 1:nrow(org.gMicrobData)) {
    # extract positions of positive taxa
    taxa.nonzero.pos=which(org.gMicrobData[i,]!=0)
    # only consider subjects with >= 2 positive taxa
    if (length(taxa.nonzero.pos) >= 2) {
      twoPos.vec[i] <- TRUE
      pairs=combn(x=sort(taxa.nonzero.pos,decreasing=T),m=2)
      k2k1.lst[[i]] <- pairs
      
      # create Ai, then calculate AZi
      Ai <- matrix(0, nrow = ncol(pairs), ncol = K)
      Ai[cbind(1:ncol(pairs), pairs[2,])] <- 1 # (l,kl)
      Ai[cbind((1:ncol(pairs))[pairs[1,] < (K+1)], pairs[1,(pairs[1,] < (K+1))])] <- -1
      AZ.lst[[i]] <- t(kronecker(t(Ai), t(CovDataWithIntcp[i, , drop = FALSE])))
      
      # Wi
      W.lst[[i]] <- as.numeric(1/((1-exp(-eta*(org.gMicrobData[i,pairs[2,]])))*(1-exp(-eta*(org.gMicrobData[i,pairs[1,]])))))
      
      # Ui
      U.lst[[i]] <- as.numeric(log(org.gMicrobData[i,pairs[2,]] / org.gMicrobData[i,pairs[1,]]))
    }
  }
  twoPosSub <- which(twoPos.vec)
  
  # store parameters that won't change in loop
  results$k2k1.lst <- k2k1.lst
  results$AZ.lst <- AZ.lst
  results$W.lst <- W.lst
  results$U.lst <- U.lst
  results$twoPosSub <- twoPosSub
  
  # update alpha and phi
  ## update alpha using the given initial value of phi
  for (i in twoPosSub) {
    pairs <- k2k1.lst[[i]]
    phi_i <- phiMat.init[cbind(pairs[2,], pairs[1,])]
    AZi <- as.matrix(AZ.lst[[i]])
    Wi <- W.lst[[i]]
    Ui <- U.lst[[i]]
    
    ## fix one pair dimension issue
    if (length(Wi) == 1) {
      ZZadd <- crossprod(AZi, (Wi / phi_i) * AZi)
      ZUadd <- t(AZi)*(Wi / phi_i) * Ui
    } else {
      Wphi_diag <- diag(Wi/phi_i)
      ZZadd <- t(AZi) %*% Wphi_diag %*% AZi
      ZUadd <- t(AZi) %*% Wphi_diag %*% Ui
    }
    
    if (i==twoPosSub[1]){
      ZZsum <- ZZadd
      ZUsum <- ZUadd
    } else {
      ZZsum <- ZZsum + ZZadd
      ZUsum <- ZUsum + ZUadd
    }
  }
  # if(K == 1){
  #   ZUsum = t(ZUsum)
  # }
  # update alpha
  phi_update.mat <- matrix(as.vector(t(phiMat.init)), nrow = 1)
  if (QRdecomp) { # then use QR decomposition to solve alpha
    alpha.new <- alpha_QRdecompose(ZZsum, ZUsum)
  } else { # directly use OLS formula which contains inv(ZZsum) to solve alpha
    alpha.new <- ginv(as.matrix(ZZsum), tol = 1e-40)%*%ZUsum
  }
  if(K == 1){
    alpha.new = as.vector(alpha.new)
  }
  alpha_update.mat <- matrix(alpha.new, nrow = 1) # to store all updates of alpha
  alpha.old <- rep(-100, ncol(CovDataWithIntcp)*K)
  
  # # store parameters which are updated in loop
  # results$phiMat <- phiMat.init
  # results$ZZsum <- ZZsum # depends on phi
  # results$ZUsum <- ZUsum # depends on phi
  # results$alpha.hat <- alpha.new # depends on phi
  
  # save for jackknife steps
  if(jackknife){
    AZWAZ.lst = list()
    AZWU.lst = list()
  }
  # start updating loop
  ite = 0
  while (mean(abs(alpha.new - alpha.old)) > 1e-3 & ite < 10) {
    ite <- ite + 1
    alpha.old <- alpha.new
    
    # update phi
    ## initialize phi's numerator and denominator to zeros for this new loop of update
    phi_nume <- matrix(0, nrow = K+1, ncol = K+1)
    phi_deno <- matrix(0.01, nrow = K+1, ncol = K+1)
    phi_deno[upper.tri(phi_deno)] <- 0
    pairs_all_count <- matrix(0, nrow = K+1, ncol = K+1)
    for (i in twoPosSub) {
      pairs <- k2k1.lst[[i]]
      AZi <- AZ.lst[[i]]
      if(K == 1){
        AZi = t(AZi)
      }
      AZalpha <- AZi%*%alpha.old
      Wi <- W.lst[[i]]
      Ui <- U.lst[[i]]
      # update phi's numerator and denominator
      weighted.squares <- Wi*(Ui - AZalpha)^2
      pairs_all_count[cbind(pairs[2,], pairs[1,])] <- pairs_all_count[cbind(pairs[2,], pairs[1,])] + 1
      phi_nume[cbind(pairs[2,], pairs[1,])] <- phi_nume[cbind(pairs[2,], pairs[1,])] + weighted.squares
      phi_deno[cbind(pairs[2,], pairs[1,])] <- phi_deno[cbind(pairs[2,], pairs[1,])] + Wi
    }
    phiMat <- phi_nume / phi_deno
    
    # If the number of sample is insufficient for estimating pair, set phi to be NA
    phiMat[pairs_all_count <= (ncol(CovDataWithIntcp) + phi_control)] = NA
    # correction for phi too low
    phiMat[phiMat < 1e-4] = NA
    
    # replace all NA with mean(upper.tri)
    phiMat[is.na(phiMat)] = mean(phiMat[upper.tri(phiMat)], na.rm = TRUE)
    # phiMat = phiMat.init
    
    # update alpha using OLS
    for (i in twoPosSub) {
      pairs <- k2k1.lst[[i]]
      phi_i <- phiMat[cbind(pairs[2,], pairs[1,])]
      AZi <- AZ.lst[[i]]
      Wi <- W.lst[[i]]
      Ui <- U.lst[[i]]
      
      # correct Inf to finite value
      WoverPhi <- InfCorrect(Wi/phi_i)
      
      if (length(WoverPhi) == 1) {
        # scalar case
        AZW <- InfCorrect(t(AZi) * WoverPhi)    # t(AZi) is 60x1 for example
        AZWAZ <- InfCorrect(AZW %*% AZi)        # matrix multiply works
        AZWU <- InfCorrect(AZW * Ui)
      } else {
        # general case
        AZW <- InfCorrect(t(AZi) %*% diag(WoverPhi))
        AZWAZ <- InfCorrect(AZW %*% AZi)
        AZWU <- InfCorrect(AZW %*% Ui)
      }
      
      if(jackknife){
        AZWAZ.lst[[i]] = AZWAZ
        AZWU.lst[[i]] = AZWU
      }
      
      if (i==twoPosSub[1]){
        ZZsum <- AZWAZ
        ZUsum <- AZWU
      } else {
        ZZsum <- InfCorrect(ZZsum + AZWAZ)
        ZUsum <- InfCorrect(ZUsum + AZWU)
      }
    }
    
    
    if (QRdecomp) { # use QR decomposition
      alpha.temp <- alpha_QRdecompose(ZZsum, ZUsum) # alpha_QRdecompose will return NA is
      if (any(is.na(alpha.temp))) break # if R in QR decomposition is singular
      if (sum(abs(alpha.temp))/sum(abs(alpha.old)) > abnormal_QR_threshold) break # if alpha is diverging and getting larger
      alpha.new <- alpha.temp
    } else { # use OLS formula
      # if the pseudo-inverse cannot be properly calculated: break this while-loop and use the previous alpha.hat as final result
      ZZsum.inv <- ginv(as.matrix(ZZsum), tol = 1e-40)
      sum_diag_identity <- sum(abs(diag(ZZsum.inv%*%ZZsum)))
      if (abs(sum_diag_identity - nrow(ZZsum))/nrow(ZZsum) > abnormal_ginv_threshold) break
      # else: there is nothing wrong with ginv(ZZsum), update alpha
      alpha.new <- ZZsum.inv%*%ZUsum
    }
    if(K == 1){
      alpha.new = t(alpha.new)
    }
    phi_update.mat <- rbind(phi_update.mat, as.vector(t(phiMat)))
    alpha_update.mat <- rbind(alpha_update.mat, as.numeric(alpha.new))
    
    
    # store the latest parameter values
    results$phiMat <- phiMat
    results$ZZsum <- ZZsum # depends on phi
    results$ZUsum <- ZUsum # depends on phi
    results$alpha.hat <- as.numeric(alpha.new) # depends on phi
    if(jackknife){
      results$AZWAZ.lst <- AZWAZ.lst
      results$AZWU.lst <- AZWU.lst
    }
    
    if (printloss) {
      stop("printloss = TRUE is not supported in this package.")
    }
  }
  
  
  # store all updates of alpha and phi for later debugging / check convergence performance
  results$alpha_update_mat <- alpha_update.mat
  results$phi_update_mat <- phi_update.mat
  
  
  return(results)
}

#' LASSO estimate of the coefficients for reference taxon screening
#'
#' Builds the weighted pairwise log-ratio design for one reference taxon
#' (weights \eqn{\sqrt{W}} with unit dispersion), standardizes the design
#' columns, and fits a LASSO path with [glmnet::glmnet()]. Only the
#' coefficients of the covariate of interest are penalized. The solution at
#' the end of the path (limited by `pmax`) is returned. The taxa whose
#' coefficients are non-zero are considered associated with the covariate,
#' which makes them unsuitable as reference taxa.
#'
#' @inheritParams alpha_phi_Update_iloop_OLS
#' @param testCovInd Column index (in `CovDataWithIntcp`) of the covariate of
#'   interest.
#' @param phiMat.init,abnormal_ginv_threshold,abnormal_QR_threshold,phi_control
#'   Unused; kept for interface compatibility with
#'   [alpha_phi_Update_iloop_OLS()].
#'
#' @return A list with the coefficient vector `alpha.hat` and `alpha_update_mat`,
#'   a one-row matrix holding the same coefficients.
#' @keywords internal
alpha_phi_Update_iloop_LASSO <- function(
    gMicrobData, # N by (K+1) taxa matrix
    CovDataWithIntcp, # N by (1+Q+S) covariate matrix which includes intercept
    eta, # parameter of zero probability mechanism
    RefTaxonName, 
    testCovInd,
    phiMat.init = matrix(1, nrow = ncol(gMicrobData), ncol = ncol(gMicrobData)),
    printloss = FALSE,
    abnormal_ginv_threshold = 0.025,
    abnormal_QR_threshold = 10,
    phi_control = 2
){
  results <- list()
  force(phiMat.init)  # Ensures evaluation before parallel execution
  
  # calculate K: K+1 is the number of taxa
  K <- ncol(gMicrobData) - 1
  
  # reorganize the microbiome matrix by moving the reference taxon to the last column
  refTaxonInd <- which(colnames(gMicrobData)==RefTaxonName)
  org.gMicrobData <- cbind(gMicrobData[,-refTaxonInd], gMicrobData[,refTaxonInd])
  rm(gMicrobData)
  
  k2k1.lst <- list()
  AZ.lst <- list()
  W.lst <- list()
  U.lst <- list()
  Zadd <- list()
  Uadd <- list()
  twoPos.vec <- rep(FALSE, nrow(org.gMicrobData))
  for (i in 1:nrow(org.gMicrobData)) {
    # extract positions of positive taxa
    taxa.nonzero.pos=which(org.gMicrobData[i,]!=0)
    # only consider subjects with >= 2 positive taxa
    if (length(taxa.nonzero.pos) >= 2) {
      twoPos.vec[i] <- TRUE
      pairs=combn(x=sort(taxa.nonzero.pos,decreasing=T),m=2)
      k2k1.lst[[i]] <- pairs
      
      # create Ai, then calculate AZi
      Ai <- matrix(0, nrow = ncol(pairs), ncol = K)
      Ai[cbind(1:ncol(pairs), pairs[2,])] <- 1 # (l,kl)
      Ai[cbind((1:ncol(pairs))[pairs[1,] < (K+1)], pairs[1,(pairs[1,] < (K+1))])] <- -1
      AZ.lst[[i]] <- t(kronecker(t(Ai), t(CovDataWithIntcp[i, , drop = FALSE])))
      
      # Wi
      W.lst[[i]] <- as.numeric(1/((1-exp(-eta*(org.gMicrobData[i,pairs[2,]])))*(1-exp(-eta*(org.gMicrobData[i,pairs[1,]])))))
      
      # Ui
      U.lst[[i]] <- as.numeric(log(org.gMicrobData[i,pairs[2,]] / org.gMicrobData[i,pairs[1,]]))
    }
  }
  twoPosSub <- which(twoPos.vec)
  
  # store parameters that won't change in loop
  # results$k2k1.lst <- k2k1.lst
  # results$AZ.lst <- AZ.lst
  # results$W.lst <- W.lst
  # results$U.lst <- U.lst
  # results$twoPosSub <- twoPosSub
  
  # update alpha and phi
  ## update alpha using the given initial value of phi
  for (i in twoPosSub) {
    pairs <- k2k1.lst[[i]]
    # phi_i <- phiMat.init[cbind(pairs[2,], pairs[1,])]
    phi_i <- 1
    AZi <- as.matrix(AZ.lst[[i]])
    Wi <- W.lst[[i]]
    Ui <- U.lst[[i]]
    
    ## fix one pair dimension issue
    if (length(Wi) == 1) {
      Zadd <- sqrt(Wi / phi_i) * AZi
      Uadd <- sqrt(Wi / phi_i) * Ui
    } else {
      Wphi_diag <- sqrt(diag(Wi/phi_i))
      Zadd <- Wphi_diag %*% AZi
      Uadd <- Wphi_diag %*% Ui
    }
    
    if (i==twoPosSub[1]){
      Z_tilda <- Zadd
      U_tilda <- Uadd
    } else {
      Z_tilda <- rbind(Z_tilda, Zadd)
      U_tilda <- c(U_tilda, Uadd)
    }
  }
  # if(K == 1){
  #   ZUsum = t(ZUsum)
  # }
  # update alpha with LASSO
  penalty.factor = rep(1:ncol(CovDataWithIntcp) %in% testCovInd, K)
  Z_tilda_scaled = apply(Z_tilda, 2, function(col) {
    s = sd(col)
    if (s == 0) {
      return(col)  # leave unchanged
    } else {
      return((col - mean(col)) / s)
    }
  })
  U_tilda_scaled = scale(U_tilda)
  glm_lasso <- glmnet(Z_tilda_scaled, U_tilda, alpha = 1, 
                      penalty.factor = penalty.factor, 
                      standardize = FALSE,
                      pmax = sum(penalty.factor)*2/3 + sum(!penalty.factor))
  alpha.new <- glm_lasso$beta[, ncol(glm_lasso$beta)]
  
  if(K == 1){
    alpha.new = as.vector(alpha.new)
  }
  alpha_update.mat <- matrix(alpha.new, nrow = 1) # to store all updates of alpha
  # alpha.old <- rep(-100, ncol(CovDataWithIntcp)*K)
  
  # # store parameters which are updated in loop
  # results$phiMat <- phiMat.init
  # results$ZZsum <- ZZsum # depends on phi
  # results$ZUsum <- ZUsum # depends on phi
  results$alpha.hat <- alpha.new # depends on phi
  # store all updates of alpha and phi for later debugging / check convergence performance
  results$alpha_update_mat <- alpha_update.mat
  # results$phi_update_mat <- phi_update.mat
  
  
  return(results)
}
