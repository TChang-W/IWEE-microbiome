#' Cross-validation job for eta selection
#'
#' Worker function run in parallel by [IWEE()] during eta selection. For one
#' reference taxon and one candidate \eqn{\eta}, it fits the model on the
#' training folds with [alpha_phi_Update_iloop_OLS()] and returns the
#' phi-weighted mean squared prediction error on the validation fold.
#'
#' @param grouptaxaInd.lst Unused; kept for interface compatibility with the
#'   other job functions.
#' @param RefInd Integer index (column of `MicrobData`) of the reference
#'   taxon.
#' @param grouprefInd.lst Unused; kept for interface compatibility.
#' @param etaInd Index into `etaCand` of the \eqn{\eta} to evaluate.
#' @param etaCand Numeric vector of candidate \eqn{\eta} values.
#' @param groupfoldInd.lst List of integer vectors with the row indices of
#'   each cross-validation fold.
#' @param val_fold Index of the fold used for validation.
#' @param MicrobData Numeric matrix of taxa abundances (samples in rows)
#'   whose columns are the taxa in the analysis, with cleaned names.
#' @param CovData Numeric matrix of covariates without intercept, with
#'   cleaned names (`x1`, `x2`, ...).
#' @param testCovNewName Cleaned name of the covariate of interest.
#' @param Newton_Raphson Logical. Must be `FALSE`; the Newton-Raphson
#'   updater is not part of the package.
#' @param QRdecomp Logical; solve with a QR decomposition (`TRUE`) or
#'   [MASS::ginv()] (`FALSE`).
#'
#' @return A list with the evaluated `eta`, the validation fold
#'   (`test_fold`), the validation error `MSE_eta`, the taxa names
#'   (`group_taxa_names`), the reference taxon (`ref_taxon_name`), the
#'   coefficient vector `alpha_OLS`, the dispersion matrix `phi_hat`, and the
#'   coefficients of the covariate of interest `beta_OLS`.
#' @keywords internal
parallel_job_cv <- function(
    grouptaxaInd.lst, # stores taxa indices 
    RefInd, # which reference in this group
    grouprefInd.lst, # stores reference taxa 
    etaInd, # which eta candidate
    etaCand, # all eta candidates (vector)
    groupfoldInd.lst, # stores CV folds
    val_fold, # test fold index
    MicrobData, # the whole taxa data (matrix)
    CovData, # covariates matrix without intercept
    testCovNewName="x1",
    Newton_Raphson=FALSE,
    QRdecomp=TRUE
){
  
  # parallel computing: small jobs: (taxa group g, ref taxon r_g, eta candidate)
  gMicrobData <- MicrobData
  
  eta <- etaCand[etaInd]
  
  # prepare dataset for later analysis
  CovDataWithIntcp <- cbind(1, CovData)
  testCovInd <- 1 + which(colnames(CovData)==testCovNewName)
  RefTaxonName <- colnames(MicrobData)[RefInd]
  
  test_index = groupfoldInd.lst[[val_fold]]
  # validation
  gMicrobData_val = gMicrobData[test_index, ]
  CovDataWithIntcp_val = CovDataWithIntcp[test_index, ]
  # training
  gMicrobData_train = gMicrobData[-test_index, ]
  CovDataWithIntcp_train = CovDataWithIntcp[-test_index, ]
  
  if (Newton_Raphson) {
    stop("Newton_Raphson = TRUE is not supported in this package.")
  } else {
    # OLS update
    alpha_phi_results <- alpha_phi_Update_iloop_OLS(gMicrobData = gMicrobData_train, 
                                                    CovDataWithIntcp = CovDataWithIntcp_train, 
                                                    eta = eta,
                                                    RefTaxonName = RefTaxonName, 
                                                    QRdecomp = QRdecomp,
                                                    printloss = FALSE)
    alpha.hat <- alpha_phi_results$alpha.hat
  }
  
  
  
  phiMat.hat <- alpha_phi_results$phiMat
  phiMat.hat[is.na(phiMat.hat)] = mean(phiMat.hat[upper.tri(phiMat.hat)], na.rm = T)
  
  k2k1.lst <- alpha_phi_results$k2k1.lst
  W.lst <- alpha_phi_results$W.lst
  AZ.lst <- alpha_phi_results$AZ.lst
  U.lst <- alpha_phi_results$U.lst
  twoPosSub <- alpha_phi_results$twoPosSub
  
  rm(alpha_phi_results)
  # coefficients of interest
  beta.mat <- matrix(alpha.hat, nrow = ncol(CovDataWithIntcp))[testCovInd,]
  
  group_taxa_names = colnames(gMicrobData)
  names(beta.mat) <- setdiff(group_taxa_names, RefTaxonName)
  
  # k-fold validation MSE
  org.gMicrobData_val <- cbind(gMicrobData_val[,-RefInd], gMicrobData_val[,RefInd])
  K <- ncol(org.gMicrobData_val) - 1
  # get the parameters in validation set
  twoPos.vec_val <- rep(FALSE, nrow(org.gMicrobData_val))
  SE_eta <- c()
  for (i in 1:nrow(org.gMicrobData_val)) {
    # extract positions of positive taxa
    taxa.nonzero.pos=which(org.gMicrobData_val[i,]!=0)
    # only consider subjects with >= 2 positive taxa
    if (length(taxa.nonzero.pos) >= 2) {
      twoPos.vec_val[i] <- TRUE
      pairs=combn(x=sort(taxa.nonzero.pos,decreasing=T),m=2)
      # create Ai, then calculate AZi
      Ai <- matrix(0, nrow = ncol(pairs), ncol = K)
      Ai[cbind(1:ncol(pairs), pairs[2,])] <- 1 # (l,kl)
      Ai[cbind((1:ncol(pairs))[pairs[1,] < (K+1)], pairs[1,(pairs[1,] < (K+1))])] <- -1
      AZi <- t(kronecker(t(Ai), t(CovDataWithIntcp_val[i,])))
      # Ui
      Ui <- as.numeric(log(org.gMicrobData_val[i,pairs[2,]] / org.gMicrobData_val[i,pairs[1,]]))
      
      # include phi_i
      phi_i <- phiMat.hat[cbind(pairs[2,], pairs[1,])]
      SE_eta[i] <- mean((Ui - AZi%*%alpha.hat)^2/phi_i, na.rm = TRUE)
      
      # ignore phi_i
      # SE_eta[i] <- mean((Ui - AZi%*%alpha.hat)^2, na.rm = TRUE)
    }
  }
  
  return(list(eta = eta,
              test_fold = val_fold, 
              MSE_eta = mean(SE_eta, na.rm = TRUE),
              # selection = selection,
              group_taxa_names = colnames(gMicrobData),
              ref_taxon_name = RefTaxonName,
              alpha_OLS = alpha.hat,
              phi_hat = phiMat.hat,
              # alpha_phi_results = alpha_phi_results,
              beta_OLS = beta.mat
              # beta_p_val = beta.p_val
  ))
  
  
}

#' Group-wise cross-validation job for eta selection (experimental)
#'
#' Alternative to [parallel_job_cv()] that works on one taxa group and
#' computes a validation error averaged by taxon pair, with the number of
#' unobserved pairs weighted by the mechanism-zero probability `pCand`.
#' It is not used by [IWEE()].
#'
#' @param g Index of the taxa group.
#' @param grouptaxaInd.lst List with the taxa column indices of each group.
#' @param RefInd Index of the reference taxon within `grouprefInd.lst[[g]]`.
#' @param grouprefInd.lst List with the reference taxa column indices of each
#'   group.
#' @param etaInd Index into `etaCand` (and `pCand`) of the candidate to
#'   evaluate.
#' @param etaCand Numeric vector of candidate \eqn{\eta} values.
#' @param pCand Numeric vector of mechanism-zero proportions matching
#'   `etaCand`.
#' @param groupfoldInd.lst List of row indices of each cross-validation fold.
#' @param val_fold Index of the validation fold.
#' @param MicrobData Numeric matrix of taxa abundances.
#' @param CovData Numeric matrix of covariates without intercept.
#' @param testCovNewName Cleaned name of the covariate of interest.
#' @param Newton_Raphson Logical. Must be `FALSE`.
#' @param QRdecomp Logical; use a QR decomposition instead of
#'   [MASS::ginv()].
#'
#' @return A list with the same components as [parallel_job_cv()] plus the
#'   group index `group`.
#' @keywords internal
parallel_job_cv_v2 <- function(
    g, # which group
    grouptaxaInd.lst, # stores taxa indices for all groups
    RefInd, # which reference in this group
    grouprefInd.lst, # stores reference taxa for all groups
    etaInd, # which eta candidate
    etaCand, # all eta candidates (vector)
    pCand, # mechnisam zero percentage candicates
    groupfoldInd.lst, # stores CV folds for all groups
    val_fold, # test fold index
    MicrobData, # the whole taxa data (matrix)
    CovData, # covariates matrix without intercept
    testCovNewName="x1",
    Newton_Raphson=FALSE,
    QRdecomp=FALSE
){
  
  # parallel computing: small jobs: (taxa group g, ref taxon r_g, eta candidate)
  gMicrobData <- MicrobData[,grouptaxaInd.lst[[g]]]
  
  eta <- etaCand[etaInd]
  mech_prob = pCand[etaInd]
  
  # prepare dataset for later analysis
  CovDataWithIntcp <- cbind(1, CovData)
  testCovInd <- 1 + which(colnames(CovData)==testCovNewName)
  RefTaxonName <- colnames(MicrobData)[grouprefInd.lst[[g]]][RefInd]
  
  test_index = groupfoldInd.lst[[val_fold]]
  # validation
  gMicrobData_val = gMicrobData[test_index, ]
  CovDataWithIntcp_val = CovDataWithIntcp[test_index, ]
  # training
  gMicrobData_train = gMicrobData[-test_index, ]
  CovDataWithIntcp_train = CovDataWithIntcp[-test_index, ]
  
  if (Newton_Raphson) {
    stop("Newton_Raphson = TRUE is not supported in this package.")
  } else {
    # OLS update
    alpha_phi_results <- alpha_phi_Update_iloop_OLS(gMicrobData = gMicrobData_train, CovDataWithIntcp = CovDataWithIntcp_train, eta = eta,
                                                    RefTaxonName = RefTaxonName, 
                                                    QRdecomp = QRdecomp,
                                                    printloss = FALSE)
    alpha.hat <- alpha_phi_results$alpha.hat
  }
  
  
  
  phiMat.hat <- alpha_phi_results$phiMat
  phiMat.hat[is.na(phiMat.hat)] = mean(phiMat.hat[upper.tri(phiMat.hat)], na.rm = T)
  
  k2k1.lst <- alpha_phi_results$k2k1.lst
  W.lst <- alpha_phi_results$W.lst
  AZ.lst <- alpha_phi_results$AZ.lst
  U.lst <- alpha_phi_results$U.lst
  twoPosSub <- alpha_phi_results$twoPosSub
  
  
  # coefficients of interest
  beta.mat <- matrix(alpha.hat, nrow = ncol(CovDataWithIntcp))[testCovInd,]
  
  group_taxa_names = colnames(gMicrobData)
  names(beta.mat) <- setdiff(group_taxa_names, RefTaxonName)
  
  # k-fold validation MSE
  org.gMicrobData_val <- cbind(gMicrobData_val[,-RefInd], gMicrobData_val[,RefInd])
  K <- ncol(org.gMicrobData_val) - 1
  
  pairs_count = matrix(0, K+1, K+1)
  pairs_SE_sum = matrix(0, K+1, K+1)
  for (i in 1:nrow(org.gMicrobData_val)) {
    # extract positions of positive taxa
    taxa.nonzero.pos=which(org.gMicrobData_val[i,]!=0)
    # only consider subjects with >= 2 positive taxa
    if (length(taxa.nonzero.pos) >= 2) {
      # twoPos.vec_val[i] <- TRUE
      pairs=combn(x=sort(taxa.nonzero.pos,decreasing=T),m=2)
      # create Ai, then calculate AZi
      Ai <- matrix(0, nrow = ncol(pairs), ncol = K)
      Ai[cbind(1:ncol(pairs), pairs[2,])] <- 1 # (l,kl)
      Ai[cbind((1:ncol(pairs))[pairs[1,] < (K+1)], pairs[1,(pairs[1,] < (K+1))])] <- -1
      AZi <- t(kronecker(t(Ai), t(CovDataWithIntcp_val[i,])))
      # Ui
      Ui <- as.numeric(log(org.gMicrobData_val[i,pairs[2,]] / org.gMicrobData_val[i,pairs[1,]]))
      # count of pairs + 1
      pairs_count[cbind(pairs[2,], pairs[1,])] = pairs_count[cbind(pairs[2,], pairs[1,])] + 1
      # squared error
      pairs_SE_sum[cbind(pairs[2,], pairs[1,])] = pairs_SE_sum[cbind(pairs[2,], pairs[1,])] + (Ui - AZi%*%alpha.hat)^2
    }
  }
  pairs_ave_SE = pairs_SE_sum/(pairs_count + (nrow(org.gMicrobData_val) - pairs_count)*mech_prob)/phiMat.hat
  
  return(list(eta = eta,
              group = g,
              test_fold = val_fold, 
              MSE_eta = mean(pairs_ave_SE, na.rm = TRUE),
              # selection = selection,
              group_taxa_names = colnames(gMicrobData),
              ref_taxon_name = RefTaxonName,
              alpha_OLS = alpha.hat,
              phi_hat = phiMat.hat,
              # alpha_phi_results = alpha_phi_results,
              beta_OLS = beta.mat
              # beta_p_val = beta.p_val
  ))
  
  
}

#' LASSO job for reference taxon selection
#'
#' Worker function run in parallel by [IWEE()] in the reference-selection
#' phase. For one taxa group and one candidate reference taxon it fits a
#' penalized (LASSO) pairwise log-ratio regression with
#' [alpha_phi_Update_iloop_LASSO()].
#'
#' @param g Index of the taxa group.
#' @param grouptaxaInd.lst List with the taxa column indices of each group.
#' @param RefInd Index of the candidate reference taxon within
#'   `grouprefInd.lst[[g]]`.
#' @param grouprefInd.lst List with the candidate reference taxa column
#'   indices of each group.
#' @param etaInd Index into `etaCand` of the selected \eqn{\eta}.
#' @param etaCand Numeric vector of candidate \eqn{\eta} values.
#' @param groupfoldInd.lst Unused; kept for interface compatibility.
#' @param MicrobData Numeric matrix of taxa abundances.
#' @param CovData Numeric matrix of covariates without intercept.
#' @param testCovNewName Cleaned name of the covariate of interest.
#' @param Newton_Raphson Logical. Must be `FALSE`.
#'
#' @return A list with `eta`, `group`, the taxa names of the group
#'   (`group_taxa_names`), the reference taxon (`ref_taxon_name`), the full
#'   coefficient vector `alpha_OLS`, and the coefficients of the covariate of
#'   interest `beta_OLS`.
#' @keywords internal
parallel_job_ref <- function(
    g, # which group
    grouptaxaInd.lst, # stores taxa indices for all groups
    RefInd, # which reference in this group
    grouprefInd.lst, # stores reference taxa for all groups
    etaInd, # which eta candidate
    etaCand, # all eta candidates (vector)
    groupfoldInd.lst, # stores CV folds for all groups
    MicrobData, # the whole taxa data (matrix)
    CovData, # covariates matrix without intercept
    testCovNewName="x1",
    Newton_Raphson=FALSE
){
  # parallel computing: small jobs: (taxa group g, ref taxon r_g, eta candidate)
  gMicrobData <- MicrobData[,grouptaxaInd.lst[[g]]]
  
  eta <- etaCand[etaInd]
  
  # prepare dataset for later analysis
  CovDataWithIntcp <- cbind(1, CovData)
  testCovInd <- 1 + which(colnames(CovData)==testCovNewName)
  RefTaxonName <- colnames(MicrobData)[grouprefInd.lst[[g]]][RefInd]
  
  if (Newton_Raphson) {
    stop("Newton_Raphson = TRUE is not supported in this package.")
  } else {
    # LASSO update
    alpha_phi_results <- alpha_phi_Update_iloop_LASSO(gMicrobData = gMicrobData, CovDataWithIntcp = CovDataWithIntcp, eta = eta,
                                                      RefTaxonName = RefTaxonName, testCovInd = testCovInd,
                                                      printloss = FALSE)
    alpha.hat <- alpha_phi_results$alpha.hat
  }
  
  
  
  # phiMat.hat <- alpha_phi_results$phiMat
  # phiMat.hat[is.na(phiMat.hat)] = mean(phiMat.hat[upper.tri(phiMat.hat)], na.rm = T)
  
  # coefficients of interest
  beta.mat <- matrix(alpha.hat, nrow = ncol(CovDataWithIntcp))[testCovInd,]
  
  
  return(list(eta = eta,
              group = g,
              # selection = selection,
              group_taxa_names = colnames(gMicrobData),
              ref_taxon_name = RefTaxonName,
              alpha_OLS = alpha.hat,
              # phi_hat = phiMat.hat,
              # alpha_JK = alpha.JK,
              # var_alpha_JK = var.JK,
              beta_OLS = beta.mat
              # beta_JK = beta.JK.mat,
              # beta_var = var.beta.mat,
              # beta_p_val = beta.p_val
  ))
  
  
}


#' Jackknife estimation job for final coefficient estimates
#'
#' Worker function run in parallel by [IWEE()] in the final phase. For one
#' taxa group and one final reference taxon it fits the model on all samples
#' with [alpha_phi_Update_iloop_OLS()] and then applies a block jackknife to
#' obtain bias-corrected coefficients and their variances.
#'
#' @param g Index of the taxa group.
#' @param grouptaxaInd.lst List with the taxa column indices of each group.
#' @param RefInd Index of the reference taxon within `grouprefInd.lst[[g]]`.
#' @param grouprefInd.lst List with the final reference taxa column indices
#'   of each group.
#' @param etaInd Index into `etaCand` of the selected \eqn{\eta}.
#' @param etaCand Numeric vector of candidate \eqn{\eta} values.
#' @param MicrobData Numeric matrix of taxa abundances.
#' @param CovData Numeric matrix of covariates without intercept.
#' @param testCovNewName Cleaned name of the covariate of interest.
#' @param Newton_Raphson Logical. Must be `FALSE`.
#' @param QRdecomp Logical; use a QR decomposition instead of
#'   [MASS::ginv()].
#' @param jackknife_block Number of jackknife blocks. Must not exceed the
#'   number of samples with at least two positive taxa.
#'
#' @return A list with `eta`, `group`, `group_taxa_names`, `ref_taxon_name`,
#'   the jackknife coefficients of the covariate of interest (`beta_JK`),
#'   their variances (`beta_var`, zeros replaced by `1e-20`), z-test p-values
#'   (`beta_p_val`), and the underlying fit (`alpha_phi_results`).
#' @keywords internal
parallel_job_jk <- function(
    g, # which group
    grouptaxaInd.lst, # stores taxa indices for all groups
    RefInd, # which reference in this group
    grouprefInd.lst, # stores reference taxa for all groups
    etaInd, # which eta candidate
    etaCand, # all eta candidates (vector)
    MicrobData, # the whole taxa data (matrix)
    CovData, # covariates matrix without intercept
    testCovNewName="x1",
    Newton_Raphson=FALSE,
    QRdecomp=FALSE,
    jackknife_block = jackknife_block
){
  
  # parallel computing: small jobs: (taxa group g, ref taxon r_g, eta candidate)
  gMicrobData <- MicrobData[,grouptaxaInd.lst[[g]]]
  
  eta <- etaCand[etaInd]
  
  # prepare dataset for later analysis
  CovDataWithIntcp <- cbind(1, CovData)
  testCovInd <- 1 + which(colnames(CovData)==testCovNewName)
  RefTaxonName <- colnames(MicrobData)[grouprefInd.lst[[g]]][RefInd]
  
  if (Newton_Raphson) {
    stop("Newton_Raphson = TRUE is not supported in this package.")
  } else {
    # OLS update
    alpha_phi_results <- alpha_phi_Update_iloop_OLS(gMicrobData = gMicrobData, CovDataWithIntcp = CovDataWithIntcp, eta = eta,
                                                    RefTaxonName = RefTaxonName, 
                                                    QRdecomp = QRdecomp,
                                                    printloss = FALSE,
                                                    jackknife = TRUE)
    if(is.null(alpha_phi_results$alpha.hat)){
      alpha.hat <- alpha_phi_results$alpha_update_mat
    } else {
      alpha.hat <- alpha_phi_results$alpha.hat
    }
    
    ZZadd <- alpha_phi_results$AZWAZ.lst
    ZUadd <- alpha_phi_results$AZWU.lst
  }
  
  
  
  phiMat.hat <- alpha_phi_results$phiMat
  k2k1.lst <- alpha_phi_results$k2k1.lst
  W.lst <- alpha_phi_results$W.lst
  AZ.lst <- alpha_phi_results$AZ.lst
  U.lst <- alpha_phi_results$U.lst
  twoPosSub <- alpha_phi_results$twoPosSub
  
  ZZsum <- alpha_phi_results$ZZsum
  ZUsum <- alpha_phi_results$ZUsum
  
  # coefficients of interest
  beta.mat <- matrix(alpha.hat, nrow = ncol(CovDataWithIntcp))[testCovInd,]
  
  
  if(jackknife_block > length(twoPosSub)) stop("Jackknife block size too big.")
  
  alpha.JK.mat <- matrix(0, nrow = jackknife_block, ncol = length(alpha.hat))
  block_size =  rep(floor(length(twoPosSub)/jackknife_block), jackknife_block)
  
  remain = length(twoPosSub) %% jackknife_block
  if (remain > 0) {
    block_size[1:remain] = block_size[1:remain] + 1
  }
  blocks = split(sample(twoPosSub), rep(1:jackknife_block, times = block_size))
  for (jk_i in 1:length(blocks)) {
    ZZ_i = ZZsum
    ZU_i = ZUsum
    for (i in blocks[[jk_i]]) {
      # subInd <- twoPosSub[i]
      ZZ_i <- ZZ_i - ZZadd[[i]]
      ZU_i <- ZU_i - ZUadd[[i]]
    }
    
    if (QRdecomp) {
      alpha.JK.mat[jk_i,] <- as.numeric(alpha_QRdecompose(ZZ_i, ZU_i))
    } else {
      alpha.JK.mat[jk_i,] <- ginv(as.matrix(ZZ_i), tol = 1e-50)%*%ZU_i
    }
    
  }
  
  ### JK alpha estimate
  alpha.JK <- alpha.hat - (jackknife_block-1)*colMeans(alpha.JK.mat - matrix(rep(alpha.hat, jackknife_block), nrow = jackknife_block, byrow = TRUE))
  beta.JK.mat <- matrix(alpha.JK, nrow = ncol(CovDataWithIntcp))[testCovInd,]
  ### Var_JK for p-value
  alpha.JK.mean <- colMeans(alpha.JK.mat)
  alpha.JK.diff = sweep(alpha.JK.mat, 2, alpha.JK.mean, FUN = "-")
  # alpha.JK.diff <- alpha.JK.mat - matrix(rep(alpha.JK, jackknife_block), nrow = jackknife_block, byrow = TRUE)
  var.JK <- (jackknife_block-1)/jackknife_block * t(alpha.JK.diff)%*%alpha.JK.diff
  
  ## select taxa dependent on covariates of interest: Zg
  var.beta.mat <- matrix(diag(var.JK), nrow = ncol(CovDataWithIntcp))[testCovInd,]
  var.beta.mat[var.beta.mat==0] <- 1e-20 # since 0 was observed in variance and caused NaN issue for setting1 scenario1 with seed77
  
  ## z-test: p-value
  beta.p_val <- pnorm(abs(beta.JK.mat/sqrt(var.beta.mat)), lower.tail = FALSE)*2
  # selection <- as.numeric(beta.p_val < quantile(beta.p_val, probs=1/3))
  
  
  return(list(eta = eta,
              group = g,
              # MSE_eta = mean(SE_eta),
              # selection = selection,
              group_taxa_names = colnames(gMicrobData),
              ref_taxon_name = RefTaxonName,
              # alpha_OLS = alpha.hat,
              # phi_hat = phiMat.hat,
              # alpha_JK = alpha.JK,
              # var_alpha_JK = var.JK,
              # beta_OLS = beta.mat,
              beta_JK = beta.JK.mat,
              beta_var = var.beta.mat,
              beta_p_val = beta.p_val,
              alpha_phi_results = alpha_phi_results
  ))
  
  
}
