#' Inverse-weighted estimating equation (IWEE) analysis of microbiome data
#'
#' Estimates the association between microbial taxa abundances and covariates
#' of interest from zero-inflated, compositional microbiome data using the
#' inverse-weighted estimating equation (IWEE) approach. Zeros are treated as
#' arising from a mixture of true (biological) and mechanism (sampling-depth)
#' zeros; the zero mechanism is governed by a parameter \eqn{\eta}.
#'
#' The procedure has four phases:
#' \enumerate{
#'   \item **Eta selection.** A grid of `nEta` candidate values of \eqn{\eta},
#'     corresponding to mechanism-zero proportions from 1 down to 0, is
#'     computed with [estiEta()]. Up to 20 taxa with 40-50% sparsity (or, if
#'     fewer than 10 such taxa exist, the 10 least sparse taxa) are used as
#'     reference taxa, and the \eqn{\eta} with the smallest cross-validated
#'     prediction error (summed over the first four of these reference taxa)
#'     is chosen.
#'   \item **Group partition and reference selection.** Taxa are randomly
#'     partitioned into groups of at least `Ng` taxa. In each group, `nRef`
#'     randomly chosen candidate reference taxa are screened with a LASSO
#'     fit (see [alpha_phi_Update_iloop_LASSO()]), and the `nBestRef` taxa
#'     with the fewest non-zero LASSO coefficients across candidates (i.e.
#'     the least likely to be associated with the covariate) are kept as
#'     final reference taxa.
#'   \item **Final estimation.** For each group and final reference taxon the
#'     coefficients are estimated by iteratively reweighted least squares
#'     (see [alpha_phi_Update_iloop_OLS()]) and bias-corrected and
#'     variance-estimated by a block jackknife.
#'   \item **Averaging.** Estimates and variances are averaged over the
#'     reference taxa of each group. Wald-type z-test p-values are computed,
#'     adjusted with the Benjamini-Hochberg (FDR) method, and taxa with
#'     adjusted p-value below 0.05 are flagged as significant.
#' }
#' Phases 1 to 3 run in parallel with a PSOCK cluster created through
#' \pkg{doParallel}. Progress and timing messages are printed with [cat()].
#'
#' @param MicrobData Microbiome data: a matrix or data frame with one row per
#'   sample, one column per taxon, and an additional sample-identifier column
#'   named by `linkIDname`. Alternatively, a path to a `.csv` or `.tsv` file
#'   with this layout. Entries are non-negative abundances (counts); zeros are
#'   allowed.
#' @param CovData Covariate data: a matrix or data frame (or `.csv`/`.tsv`
#'   file path) with one row per sample, the identifier column `linkIDname`,
#'   and one column per covariate.
#' @param num_cores Number of worker processes for the parallel stages. If
#'   `NA` (the default), the number of physical cores reported by
#'   [parallel::detectCores()] is used.
#' @param testCov Character vector with the name(s) of the covariate(s) of
#'   interest, as they appear in the column names of `CovData`. Note that the
#'   estimation stages use a single test covariate.
#' @param ctrlCov Character vector with the names of confounders to adjust
#'   for. Defaults to `paste0("W", 1:2)`. If it is empty and `ctrlMany` is
#'   `TRUE`, all covariates other than `testCov` are used as controls.
#' @param linkIDname Name of the sample-identifier column present in both
#'   `MicrobData` and `CovData`, used to match the two tables.
#' @param QRdecomp Logical. If `TRUE`, linear systems are solved with a QR
#'   decomposition (see [alpha_QRdecompose()]); if `FALSE`, with the
#'   Moore-Penrose pseudo-inverse [MASS::ginv()]. The LASSO reference-selection
#'   step does not depend on this choice.
#' @param nGroup Ignored. The number of groups is always recomputed as
#'   `ncol(MicrobData) %/% Ng`, where the column count is taken after
#'   removing taxa without reads. The argument is retained for backward
#'   compatibility.
#' @param nEta Number of candidate values of \eqn{\eta}, spanning
#'   mechanism-zero proportions from 1 to 0 in equal steps.
#' @param Ng Minimum number of taxa per group after partitioning. The last
#'   group absorbs the remaining taxa. Must not exceed the number of taxa.
#' @param nRef Number of candidate reference taxa drawn at random in each
#'   group. Must not exceed `Ng`.
#' @param nBestRef Number of final reference taxa retained per group.
#' @param k_fold Number of cross-validation folds. The current implementation
#'   evaluates only the first fold when selecting \eqn{\eta}.
#' @param jackknife_block Number of blocks used by the block jackknife. It
#'   must not exceed the number of samples with at least two positive taxa.
#' @param testMany Logical. If `TRUE` and `testCov` is empty, all covariates
#'   are treated as covariates of interest.
#' @param ctrlMany Logical. If `TRUE` and `ctrlCov` is empty, all covariates
#'   other than `testCov` are treated as confounders.
#' @param randseed Integer random seed (see [set.seed()]), controlling the
#'   group partition, reference candidates, folds, and jackknife blocks.
#'
#' @return A list with components:
#' \describe{
#'   \item{`eta_cvMSE`}{List (one element per reference taxon used in the eta
#'     selection) of matrices with the cross-validation error for each
#'     candidate \eqn{\eta}, and its sum over folds.}
#'   \item{`etaMSE`}{List of the corresponding error vectors over the
#'     candidate \eqn{\eta} values.}
#'   \item{`etaCand`}{Numeric vector of the `nEta` candidate \eqn{\eta}
#'     values.}
#'   \item{`bestEtaCount`}{Integer vector: for each candidate \eqn{\eta}, the
#'     number of reference taxa for which it gave the smallest error.}
#'   \item{`bestEtaInd`}{Index (into `etaCand`) of the selected \eqn{\eta}.}
#'   \item{`final_ref_names.lst`}{List with one character vector per group of
#'     the selected reference taxa (internal names `microb1`, `microb2`,
#'     ...).}
#'   \item{`final_beta_var`}{Named numeric vector of jackknife variances of
#'     the coefficient estimates, averaged over reference taxa.}
#'   \item{`final_beta_est`}{Named numeric vector of jackknife bias-corrected
#'     coefficient estimates for `testCov`, averaged over reference taxa.}
#'   \item{`final_beta_pval`}{Named numeric vector of two-sided z-test
#'     p-values.}
#'   \item{`final_beta_fdr`}{Named numeric vector of FDR-adjusted p-values.}
#'   \item{`final_sig_beta`}{Named logical vector, `TRUE` where
#'     `final_beta_fdr < 0.05`.}
#'   \item{`randseed`}{The random seed used.}
#'   \item{`phase1_count_coefMagnit`}{List, one data frame per group, with
#'     the number of candidate reference fits in which each taxon had a
#'     non-zero LASSO coefficient, ordered increasingly (the first
#'     `nBestRef` rows give the final reference taxa).}
#' }
#' The `final_beta_*` and `final_sig_beta` vectors are ordered by taxon and
#' named with the original taxon names found in `MicrobData`.
#'
#' @seealso [metaData()] for the data preparation step and [estiEta()] for
#'   the computation of the \eqn{\eta} candidates.
#'
#' @examples
#' # Example data: 100 samples, 400 taxa, binary covariate X1 and confounders
#' # W1, W2 (simulated; see ?load_iwee_example).
#' dat <- load_iwee_example()
#'
#' \donttest{
#' # A reduced run (first 100 taxa, small tuning grid) to keep run time short
#' res <- IWEE(
#'   MicrobData = dat$MicrobData[, 1:101],
#'   CovData = dat$CovData,
#'   num_cores = 2,
#'   testCov = "X1",
#'   ctrlCov = c("W1", "W2"),
#'   linkIDname = "id",
#'   QRdecomp = TRUE,
#'   nEta = 3,
#'   Ng = 50,
#'   nRef = 10,
#'   nBestRef = 2,
#'   k_fold = 2,
#'   testMany = FALSE,
#'   ctrlMany = TRUE,
#'   randseed = 1
#' )
#' head(res$final_beta_est)
#' names(which(res$final_sig_beta))
#' }
#' @export
IWEE <- function(
    MicrobData,
    CovData,
    num_cores = NA, 
    testCov = "X1", # covariates of interest
    ctrlCov = paste0("W", 1:2), # confounders
    linkIDname = "id",
    QRdecomp = FALSE,
    nGroup = 3,
    nEta = 11, # number of eta candidates corresponding to p = 0 to 1
    Ng = 16, # the smallest number of taxa in a taxa group after partitioning
    nRef = 5, # the number of selected reference candidates in each group
    nBestRef = 2, # the number of final selected best reference taxa in each group
    k_fold = 5, # number of folds 
    jackknife_block = 10, # jackknife size
    testMany = FALSE,
    ctrlMany = TRUE,
    randseed = 213
){
  set.seed(randseed)
  
  
  
  cat("Data cleaning started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  flush.console()  # Ensures it prints immediately in interactive environments
  # clean dataset
  runMeta <- metaData(MicrobData = MicrobData, CovData = CovData, linkIDname = linkIDname,
                      testCov = testCov, ctrlCov = ctrlCov, testMany = testMany, ctrlMany = ctrlMany)
  data <- runMeta$data
  Mprefix <- runMeta$Mprefix
  covsPrefix <- runMeta$covsPrefix
  ## original and new names for taxa
  MicrobDataOrigName <- runMeta$microbName
  MicrobDataNewName <- runMeta$newMicrobNames
  MicrobNameMap <- rbind(MicrobDataOrigName, MicrobDataNewName)
  ## new name for covariates of interest
  testCovNewName <- runMeta$testCovInNewNam
  ## obtain data with cleaned colnames and without id
  MicrobData <- data[, grepl(paste0("^",Mprefix), colnames(data)), drop = FALSE]
  CovData <- data[, grepl(paste0("^",covsPrefix), colnames(data)), drop = FALSE]
  rm(data)
  cat("Data cleaning ended at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  flush.console()  # Ensures it prints immediately in interactive environments
  
  cat("Group partition started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  flush.console()  # Ensures it prints immediately in interactive environments
  # partition taxa into smaller groups
  nTaxa <- ncol(MicrobData)
  randomShuf <- sample(nTaxa, nTaxa)
  nGroup <- nTaxa %/% Ng
  if (nTaxa < Ng) stop("the smallest number of taxa for all partitioned groups cannot be larger than the total number of taxa.")
  if (Ng < nRef) stop("the number of reference candidates cannot be larger than the number of taxa in a partitioned group.")
  grouptaxaInd.lst <- list()
  grouprefInd.lst <- list()
  if (nGroup==1) {
    grouptaxaInd.lst[[1]] <- c(1:nTaxa)
    # randomly select reference candidates for this group
    grouprefInd.lst[[1]] <- sort(sample(grouptaxaInd.lst[[1]], nRef))
  } else {
    for (g in 1:nGroup) {
      if (g != nGroup) {
        grouptaxaInd.lst[[g]] <- sort(randomShuf[((g-1)*Ng+1):(g*Ng)])
      } else {
        grouptaxaInd.lst[[g]] <- sort(randomShuf[((g-1)*Ng+1):nTaxa])
      }
      # randomly select reference candidates for this group
      grouprefInd.lst[[g]] <- sort(sample(grouptaxaInd.lst[[g]], nRef))
    }
  }
  
  N = nrow(MicrobData)
  if (N %/% k_fold <= 2)
    stop("The number of samples is too small for k-fold validation.")
  
  groupfoldInd.lst = list()
  foldrandomShuf = sample(N)  # shuffled indices
  fold_size = floor(N / k_fold)
  
  for (k in 1:k_fold) {
    if (k != k_fold) {
      inds = ((k - 1) * fold_size + 1):(k * fold_size)
    } else {
      inds = ((k - 1) * fold_size + 1):N  # include all remaining
    }
    groupfoldInd.lst[[k]] = sort(foldrandomShuf[inds])
  }
  
  
  cat("Group partition ended at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  flush.console()  # Ensures it prints immediately in interactive environments
  
  
  
  cat("Eta calculation started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  flush.console()  # Ensures it prints immediately in interactive environments
  # calculate eta's candidates
  etaCand <- c()
  pCand <- seq(1, 0, length.out = nEta)
  for (i in 1:nEta) {
    etaCand[i] <- estiEta(MicrobData, mechZeroPerc = pCand[i])
    cat("eta's Candidate: ", etaCand[i], "\n")
  }
  cat("Eta calculation ended at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  flush.console()  # Ensures it prints immediately in interactive environments
  
  
  ##################################first loop#####################################
  
  # Select up to 10 taxa with 40-50% sparsity for eta selection
  sparsity = colMeans(MicrobData == 0 | is.na(MicrobData))
  cols_target = names(sparsity[sparsity >= 0.4 & sparsity <= 0.5])
  if (length(cols_target) >= 10) {
    # If enough columns exist, randomly select up to 30
    selected_cols = sample(cols_target, min(20, length(cols_target)))
  } else {
    # select from top sparsity
    ord = order(sparsity)
    sparsity_sorted = sparsity[ord]
    upper = min(10, ncol(MicrobData))
    selected_cols = names(sparsity_sorted[1:upper])
  }
  MicrobData_selected = MicrobData[, selected_cols]
  # Print summary
  cat("Number of selected columns:", length(selected_cols), "\n")
  flush.console()
  
  cat("First set parallel tasks started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  flush.console()  # Ensures it prints immediately in interactive environments
  # foreach parallel computing for Phase 1: Association identification
  ## you can specify to use the simple OLS with ginv() or the QR decomposition with Matrix::qr() and qr.coef()
  if(is.na(num_cores)){
    num_cores <- parallel::detectCores(logical = FALSE)
  }
  cl <- makeCluster(num_cores)
  # cl <- makeForkCluster(num_cores) # even slower than makeCluster
  registerDoParallel(cl)
  cat("Allocated CPU cores:", length(cl), "\n")  
  
  num_workers <- getDoParWorkers()  # Check registered workers
  cat("Number of registered parallel workers:", num_workers, "\n")
  
  gc()
  # peakRAM(
  # system.time({ # 5s * (2*20*11/7)
  results <- foreach(RefInd = 1:4, .packages = c("MASS", "Matrix"), 
                     .export = c("parallel_job_cv", "alpha_phi_Update_iloop_OLS", "InfCorrect", "alpha_QRdecompose"),
                     .options.snow = list(preschedule = FALSE)) %:%  # Outer loop .export = c("MicrobData", "CovData", "parallel_job", "grouptaxaInd.lst")
    foreach(etaInd = 1:nEta) %:%
    foreach(val_fold = 1) %dopar% {  # Inner loop, running in parallel
      # each small job runs about 260s, 440 small jobs take about peakRAM--82600MiB
      tryCatch({
        parallel_job_cv(grouptaxaInd.lst = colnames(MicrobData_selected), testCovNewName = testCovNewName,
                        RefInd = RefInd, grouprefInd.lst = colnames(MicrobData_selected),
                        etaInd = etaInd, etaCand = etaCand, 
                        groupfoldInd.lst = groupfoldInd.lst, val_fold = val_fold,
                        MicrobData = MicrobData_selected, CovData = CovData,
                        QRdecomp = QRdecomp,
                        Newton_Raphson = FALSE)
      }, error = function(e) {
        list(error = TRUE,
             message = e$message,
             g = g,
             RefInd = RefInd,
             etaInd = etaInd)
      })
      # cat("group: ", g, ", RefInd: ", RefInd, ", etaInd: ", etaInd)
    }
  # save(results, file = "results_firstloop.RData")
  # })
  # )
  cat("First set parallel tasks ended at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  flush.console()  # Ensures it prints immediately in interactive environments
  
  
  
  cat("stopCluster(cl) started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  flush.console()  # Ensures it prints immediately in interactive environments
  # Properly stop the cluster after use
  # stopCluster(cl)
  parallel::stopCluster(cl)
  rm(cl)
  gc()  # Garbage collection to free memory
  cat("Available cores after stopCluster(cl):", parallel::detectCores(), "\n")
  flush.console()  # Ensures output is immediately written to .out file
  cat("stopCluster(cl) ended at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  flush.console()  # Ensures it prints immediately in interactive environments
  flush.console()  # Ensures it prints immediately in interactive environments
  
  ##################################eta selection#####################################
  # Elapsed_Time_sec Total_RAM_Used_MiB Peak_RAM_Used_MiB
  #         3224.874               0.6           1119953
  # save(results_OLS, file = "results_OLS.RData")
  # load("../sim0/results_OLS.RData")
  
  
  
  
  
  cat("Eta selection started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  flush.console()  # Ensures it prints immediately in interactive environments
  # selection analysis
  ## eta selection result
  bestEtaCount <- rep(0, nEta)
  
  ### to store the intermediate coefficient results: a list of nGroup elements, each element should be a vector of avg magnitude
  count_intMedCoefMag.lst <- list()
  eta_cvMSE <- list()
  etaMSE <- list()
  g_eta_MSE <- list()
  group_taxa_name <- results[[1]][[1]][[1]]$group_taxa_names
  
  beta_magnitude.mat <- matrix(NA, nrow = nRef, ncol = length(group_taxa_name))
  colnames(beta_magnitude.mat) <- group_taxa_name
  
  eta_cvMSE <- list()
  etaMSE <- list()
  
  for (RefInd in 1:4) {
    ref_taxon_name <- results[[RefInd]][[1]][[1]]$ref_taxon_name
    non_ref_taxon_name <- setdiff(group_taxa_name, ref_taxon_name)
    
    eta_cvMSE[[RefInd]] <- matrix(NA, nrow = nEta, ncol = 1)
    # look through all eta's to decide the best eta (with smallest MSE)
    for (etaInd in 1:nEta) {
      for (fold in 1) {
        MSE_eta_tmp = results[[RefInd]][[etaInd]][[fold]]$MSE_eta
        if(is.null(MSE_eta_tmp)){
          eta_cvMSE[[RefInd]][etaInd, fold] = NA
        }else{
          eta_cvMSE[[RefInd]][etaInd, fold] = MSE_eta_tmp
        }
      }
    }
    etaMSE[[RefInd]] = rowSums(eta_cvMSE[[RefInd]], na.rm = TRUE)
    eta_cvMSE[[RefInd]] = cbind(eta_cvMSE[[RefInd]], etaMSE[[RefInd]])
    colnames(eta_cvMSE[[RefInd]]) = c(paste0("fold", 1), "Sum")
    bestEtaInd = which.min(etaMSE[[RefInd]])
    # collect best_eta counts
    bestEtaCount[bestEtaInd] <- bestEtaCount[bestEtaInd] + 1
    
    g_eta_MSE = Reduce(`+`, etaMSE)
    
  }
  
  bestEtaInd = which.min(g_eta_MSE)
  cat("Eta selection ended at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  flush.console()  # Ensures it prints immediately in interactive environments
  
  ##################################second loop#####################################
  
  cat("Second set parallel tasks started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  flush.console()  # Ensures it prints immediately in interactive environments
  # foreach parallel computing for Phase 1: Association identification
  ## you can specify to use the simple OLS with ginv() or the QR decomposition with Matrix::qr() and qr.coef()
  if(is.na(num_cores)){
    num_cores <- parallel::detectCores(logical = FALSE)
  }
  cl <- makeCluster(num_cores)
  # cl <- makeForkCluster(num_cores) # even slower than makeCluster
  registerDoParallel(cl)
  cat("Allocated CPU cores:", length(cl), "\n")
  
  num_workers <- getDoParWorkers()  # Check registered workers
  cat("Number of registered parallel workers:", num_workers, "\n")
  
  gc()
  # peakRAM(
  # system.time({ # 5s * (2*20*11/7)
  results_ref <- foreach(g = 1:nGroup, .packages = c("MASS", "Matrix", "glmnet"), 
                         .export = c("parallel_job_ref", "alpha_phi_Update_iloop_LASSO", "InfCorrect", "alpha_QRdecompose"),
                         .options.snow = list(preschedule = FALSE)) %:%  # Outer loop .export = c("MicrobData", "CovData", "parallel_job", "grouptaxaInd.lst")
    foreach(RefInd = 1:nRef) %dopar% {  # Inner loop, running in parallel
      # each small job runs about 260s, 440 small jobs take about peakRAM--82600MiB
      tryCatch({
        parallel_job_ref(g = g, grouptaxaInd.lst = grouptaxaInd.lst, testCovNewName = testCovNewName,
                         RefInd = RefInd, grouprefInd.lst = grouprefInd.lst,
                         etaInd = bestEtaInd, etaCand = etaCand, 
                         groupfoldInd.lst = groupfoldInd.lst,
                         MicrobData = MicrobData, CovData = CovData,
                         Newton_Raphson = FALSE)
      }, error = function(e) {
        list(error = TRUE,
             message = e$message,
             g = g,
             RefInd = RefInd,
             etaInd = etaInd)
      })
      # cat("group: ", g, ", RefInd: ", RefInd, ", etaInd: ", etaInd)
    }
  # save(results_ref, file = "results_secondloop.RData")
  # })
  # )
  cat("Second set parallel tasks ended at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  flush.console()  # Ensures it prints immediately in interactive environments
  
  
  
  cat("stopCluster(cl) started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  flush.console()  # Ensures it prints immediately in interactive environments
  # Properly stop the cluster after use
  # stopCluster(cl)
  parallel::stopCluster(cl)
  rm(cl)
  gc()  # Garbage collection to free memory
  cat("Available cores after stopCluster(cl):", parallel::detectCores(), "\n")
  flush.console()  # Ensures output is immediately written to .out file
  cat("stopCluster(cl) ended at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  flush.console()  # Ensures it prints immediately in interactive environments
  flush.console()  # Ensures it prints immediately in interactive environments
  
  cat("Reference taxa selection started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  flush.console()  # Ensures it prints immediately in interactive environments
  # selection analysis
  ## list of final selected reference taxa: 2 final ref taxa for each group
  final_ref_names.lst <- list()
  final_ref_Ind.lst <- list()
  
  ### to store the intermediate coefficient results: a list of nGroup elements, each element should be a vector of avg magnitude
  count_intMedCoefMag.lst <- list()
  for (g in 1:nGroup) {
    group_taxa_name <- results_ref[[g]][[1]]$group_taxa_names
    
    beta_magnitude.mat <- matrix(NA, nrow = nRef, ncol = length(group_taxa_name))
    colnames(beta_magnitude.mat) <- group_taxa_name
    
    taxa_selection_count <- rep(0, length(group_taxa_name))
    names(taxa_selection_count) <- group_taxa_name
    
    for (RefInd in 1:nRef) {
      ref_taxon_name <- results_ref[[g]][[RefInd]]$ref_taxon_name
      non_ref_taxon_name <- setdiff(group_taxa_name, ref_taxon_name)
      # extract the beta estimate
      
      bestEta.beta = results_ref[[g]][[RefInd]]$beta_OLS
      beta_magnitude.mat[RefInd, non_ref_taxon_name] <- bestEta.beta
      selection <- bestEta.beta!=0
      
      taxa_selection_count[non_ref_taxon_name] <- taxa_selection_count[non_ref_taxon_name] + selection
    }
    # for (RefInd in 1:nRef) {
    #   ref_taxon_name <- results_ref[[g]][[RefInd]]$ref_taxon_name
    #   non_ref_taxon_name <- setdiff(group_taxa_name, ref_taxon_name)
    #   # given the best eta, extract the beta and phi estimate accordingly
    #   # phi estimate is the one on reference taxa
    #   bestEta.beta = NULL
    #   bestEta.phi = NULL
    #   
    #   betatmp = results_ref[[g]][[RefInd]]$beta_OLS
    #   phitmp = results_ref[[g]][[RefInd]]$phi_hat
    #   phitmp = head(phitmp[, ncol(phitmp)], -1)
    #   bestEta.beta = rbind(bestEta.beta, betatmp)
    #   bestEta.phi = rbind(bestEta.phi, phitmp)
    #   # Use adjusted beta estimate to determine selection
    #   bestEta.beta.adj = colMeans(bestEta.beta/bestEta.phi)
    #   beta_magnitude.mat[RefInd, non_ref_taxon_name] <- bestEta.beta.adj
    #   
    #   selection <- as.numeric(abs(bestEta.beta.adj) > quantile(abs(bestEta.beta.adj), probs=1/3))
    #   
    #   taxa_selection_count[non_ref_taxon_name] <- taxa_selection_count[non_ref_taxon_name] + selection
    # }
    
    # average beta.hat within each group: the avg should be biased from truth, but seems it won't influence the final selection
    avg.beta.hat <- colMeans(abs(beta_magnitude.mat), na.rm = TRUE)
    # select two most likely reference taxa for each group
    taxa_selection_beta_data <- data.frame(
      taxa = names(taxa_selection_count),  # Taxa names (assuming they are column names)
      count = as.numeric(taxa_selection_count)  # First row has selection counts
      # beta = as.numeric(avg.beta.hat)  # First row has beta estimates
    )
    taxa_selection_beta_data <- taxa_selection_beta_data[order(taxa_selection_beta_data$count), ]
    final_ref_names.lst[[g]] <- taxa_selection_beta_data$taxa[1:nBestRef]
    final_ref_Ind.lst[[g]] <- as.numeric(gsub(Mprefix, "", final_ref_names.lst[[g]]))
    
    # to store count and intermediate beta's magnitude
    count_intMedCoefMag.lst[[g]] <- taxa_selection_beta_data
  }
  cat("Reference taxa selection ended at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  flush.console()  # Ensures it prints immediately in interactive environments
  
  ##################################RESULTS#####################################
  cat("Final set parallel tasks started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  flush.console()  # Ensures it prints immediately in interactive environments
  # final regression
  # foreach parallel computing for Phase 1: Association identification
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  final_results <- foreach(g = 1:nGroup, .packages = c("MASS", "Matrix"),
                           .export = c("parallel_job_jk", "alpha_phi_Update_iloop_OLS", "InfCorrect", "alpha_QRdecompose")) %:%  # Outer loop .export = c("MicrobData", "CovData", "parallel_job", "grouptaxaInd.lst")
    foreach(RefInd = 1:nBestRef) %dopar% {  # Inner loop, running in parallel
      # each small job runs about 260s
      parallel_job_jk(g = g, grouptaxaInd.lst = grouptaxaInd.lst, testCovNewName = testCovNewName,
                      RefInd = RefInd, grouprefInd.lst = final_ref_Ind.lst,
                      etaInd = bestEtaInd, etaCand = etaCand,
                      MicrobData = MicrobData, CovData = CovData,
                      QRdecomp = QRdecomp, Newton_Raphson = FALSE,
                      jackknife_block = jackknife_block)
    }
  # Properly stop the cluster after use
  stopCluster(cl)
  cat("Final set parallel tasks ended at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  flush.console()  # Ensures it prints immediately in interactive environments
  
  
  
  
  # save(final_results, file = "final_results_OLS.RData")
  # load("../sim0/final_results_OLS.RData")
  # save(final_results_QR_qrcoef, file = "final_results_QR_qrcoef.RData")
  # load("../sim0/final_results_QR_qrcoef.RData")
  
  
  
  cat("Final Averaging started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  flush.console()  # Ensures it prints immediately in interactive environments
  # Take average of beta.hat for each group
  final_beta_est <- list()
  final_beta_var <- list()
  for (g in 1:nGroup) {
    group_taxa_name <- final_results[[g]][[1]]$group_taxa_names
    beta_JK.mat <- matrix(NA, nrow = nBestRef, ncol = length(group_taxa_name)) # 0 instead of NA for the reference taxon
    beta_var.mat <- matrix(NA, nrow = nBestRef, ncol = length(group_taxa_name)) # 0 instead of NA for the reference taxon
    colnames(beta_JK.mat) <- group_taxa_name
    colnames(beta_var.mat) <- group_taxa_name
    for (RefInd in 1:nBestRef) {
      ref_taxon_name <- final_results[[g]][[RefInd]]$ref_taxon_name
      non_ref_taxon_name <- setdiff(group_taxa_name, ref_taxon_name)
      
      beta_JK.mat[RefInd,non_ref_taxon_name] <- final_results[[g]][[RefInd]]$beta_JK
      beta_var.mat[RefInd,non_ref_taxon_name] <- final_results[[g]][[RefInd]]$beta_var
    }
    final_beta_est[[g]] <- colMeans(beta_JK.mat, na.rm = TRUE)
    final_beta_var[[g]] <- colMeans(beta_var.mat, na.rm = TRUE)
  }
  final_beta_est <- unlist(final_beta_est)
  final_beta_est <- final_beta_est[order(as.numeric(gsub(Mprefix, "", names(final_beta_est))))]
  final_beta_var <- unlist(final_beta_var)
  final_beta_var <- final_beta_var[order(as.numeric(gsub(Mprefix, "", names(final_beta_var))))]
  final_beta_pval <- pnorm(abs(final_beta_est/sqrt(final_beta_var)), lower.tail = F)*2
  final_beta_fdr = p.adjust(final_beta_pval, method = "fdr")
  final_sig_beta = final_beta_fdr < 0.05
  # map taxa names back to their original names
  name_map_vector <- setNames(MicrobNameMap[1,], MicrobNameMap[2,])
  names(final_beta_est) <- name_map_vector[names(final_beta_est)]
  names(final_beta_var) <- name_map_vector[names(final_beta_var)]
  names(final_beta_pval) <- name_map_vector[names(final_beta_pval)]
  names(final_beta_fdr) <- name_map_vector[names(final_beta_fdr)]
  names(final_sig_beta) <- name_map_vector[names(final_sig_beta)]
  cat("Final Averaging ended at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  flush.console()  # Ensures it prints immediately in interactive environments
  
  
  return(list(
    # results = results,
    # results_ref = results_ref,
    # final_results = final_results,
    eta_cvMSE = eta_cvMSE,
    etaMSE = etaMSE,
    etaCand = etaCand,
    bestEtaCount = bestEtaCount,
    bestEtaInd = bestEtaInd,
    final_ref_names.lst = final_ref_names.lst, # the selected best reference taxa
    final_beta_var = final_beta_var, # used for association identification plot
    final_beta_est = final_beta_est, # used for bias table (the true coeff values can be stored in a .csv file)
    final_beta_pval = final_beta_pval,
    final_beta_fdr = final_beta_fdr,
    final_sig_beta = final_sig_beta,
    randseed = randseed,
    phase1_count_coefMagnit = count_intMedCoefMag.lst)) # randomseed for debugging purpose
  
}
