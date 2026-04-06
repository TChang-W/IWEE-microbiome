library(MASS) # ginv()
library(Matrix) # for QR decomposition
library(foreach)
library(doParallel)
library(peakRAM)

IWEE <- function(
    MicrobData,
    CovData,
    num_cores, 
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
    testMany = F,
    ctrlMany = T,
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
  results <- foreach(g = 1:nGroup, .packages = c("MASS", "Matrix"), 
                     .export = c("parallel_job_cv", "alpha_phi_Update_iloop_OLS", "InfCorrect", "alpha_QRdecompose"),
                     .options.snow = list(preschedule = FALSE)) %:%  # Outer loop .export = c("MicrobData", "CovData", "parallel_job", "grouptaxaInd.lst")
    foreach(RefInd = 1:nRef) %:%  # Middle loop
    foreach(etaInd = 1:nEta) %:%
    foreach(val_fold = 1:k_fold ) %dopar% {  # Inner loop, running in parallel
      # each small job runs about 260s, 440 small jobs take about peakRAM--82600MiB
      tryCatch({
        parallel_job_cv(g = g, grouptaxaInd.lst = grouptaxaInd.lst, testCovNewName = testCovNewName,
                        RefInd = RefInd, grouprefInd.lst = grouprefInd.lst,
                        etaInd = etaInd, etaCand = etaCand, 
                        groupfoldInd.lst = groupfoldInd.lst, val_fold = val_fold,
                        MicrobData = MicrobData, CovData = CovData,
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
  for (g in 1:nGroup) {
    group_taxa_name <- results[[g]][[1]][[1]][[1]]$group_taxa_names
    
    beta_magnitude.mat <- matrix(NA, nrow = nRef, ncol = length(group_taxa_name))
    colnames(beta_magnitude.mat) <- group_taxa_name
    
    eta_cvMSE[[g]] <- list()
    etaMSE[[g]] <- list()
    
    for (RefInd in 1:nRef) {
      ref_taxon_name <- results[[g]][[RefInd]][[1]][[1]]$ref_taxon_name
      non_ref_taxon_name <- setdiff(group_taxa_name, ref_taxon_name)
      
      eta_cvMSE[[g]][[RefInd]] <- matrix(NA, nrow = nEta, ncol = k_fold)
      # look through all eta's to decide the best eta (with smallest MSE)
      for (etaInd in 1:nEta) {
        for (fold in 1:k_fold) {
          MSE_eta_tmp = results[[g]][[RefInd]][[etaInd]][[fold]]$MSE_eta
          if(is.null(MSE_eta_tmp)){
            eta_cvMSE[[g]][[RefInd]][etaInd, fold] = NA
          }else{
            eta_cvMSE[[g]][[RefInd]][etaInd, fold] = MSE_eta_tmp
          }
        }
      }
      etaMSE[[g]][[RefInd]] = rowSums(eta_cvMSE[[g]][[RefInd]], na.rm = TRUE)
      eta_cvMSE[[g]][[RefInd]] = cbind(eta_cvMSE[[g]][[RefInd]], etaMSE[[g]][[RefInd]])
      colnames(eta_cvMSE[[g]][[RefInd]]) = c(paste0("fold", 1:k_fold), "Sum")
      bestEtaInd = which.min(etaMSE[[g]][[RefInd]])
      # collect best_eta counts
      bestEtaCount[bestEtaInd] <- bestEtaCount[bestEtaInd] + 1
      
      g_eta_MSE[[g]] = Reduce(`+`, etaMSE[[g]])

    }
  }
  bestEtaInd = which.min(Reduce(`+`, g_eta_MSE))
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
  results_ref <- foreach(g = 1:nGroup, .packages = c("MASS", "Matrix"), 
                     .export = c("parallel_job_ref", "alpha_phi_Update_iloop_OLS", "InfCorrect", "alpha_QRdecompose"),
                     .options.snow = list(preschedule = FALSE)) %:%  # Outer loop .export = c("MicrobData", "CovData", "parallel_job", "grouptaxaInd.lst")
    foreach(RefInd = 1:nRef) %dopar% {  # Inner loop, running in parallel
      # each small job runs about 260s, 440 small jobs take about peakRAM--82600MiB
      tryCatch({
        parallel_job_ref(g = g, grouptaxaInd.lst = grouptaxaInd.lst, testCovNewName = testCovNewName,
                        RefInd = RefInd, grouprefInd.lst = grouprefInd.lst,
                        etaInd = bestEtaInd, etaCand = etaCand, 
                        groupfoldInd.lst = groupfoldInd.lst,
                        MicrobData = MicrobData, CovData = CovData,
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
      # given the best eta, extract the beta and phi estimate accordingly
      # phi estimate is the one on reference taxa
      bestEta.beta = NULL
      bestEta.phi = NULL
      for (fold in 1:k_fold) {
        betatmp = results_ref[[g]][[RefInd]]$beta_OLS
        phitmp = results_ref[[g]][[RefInd]]$phi_hat
        phitmp = head(phitmp[, ncol(phitmp)], -1)
        bestEta.beta = rbind(bestEta.beta, betatmp)
        bestEta.phi = rbind(bestEta.phi, phitmp)
      }
      
      # Use adjusted beta estimate to determine selection
      bestEta.beta.adj = colMeans(bestEta.beta/bestEta.phi)
      beta_magnitude.mat[RefInd, non_ref_taxon_name] <- bestEta.beta.adj
      
      selection <- as.numeric(abs(bestEta.beta.adj) > quantile(abs(bestEta.beta.adj), probs=1/3))
      
      taxa_selection_count[non_ref_taxon_name] <- taxa_selection_count[non_ref_taxon_name] + selection
    }
    
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
                      jackknife_block = 50)
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
    beta_JK.mat <- matrix(0, nrow = nBestRef, ncol = length(group_taxa_name)) # 0 instead of NA for the reference taxon
    beta_var.mat <- matrix(0, nrow = nBestRef, ncol = length(group_taxa_name)) # 0 instead of NA for the reference taxon
    colnames(beta_JK.mat) <- group_taxa_name
    colnames(beta_var.mat) <- group_taxa_name
    for (RefInd in 1:nBestRef) {
      ref_taxon_name <- final_results[[g]][[RefInd]]$ref_taxon_name
      non_ref_taxon_name <- setdiff(group_taxa_name, ref_taxon_name)
      
      beta_JK.mat[RefInd,non_ref_taxon_name] <- final_results[[g]][[RefInd]]$beta_JK
      beta_var.mat[RefInd,non_ref_taxon_name] <- final_results[[g]][[RefInd]]$beta_var
    }
    final_beta_est[[g]] <- colMeans(beta_JK.mat)
    final_beta_var[[g]] <- colMeans(beta_var.mat, na.rm = TRUE)
  }
  final_beta_est <- unlist(final_beta_est)
  final_beta_est <- final_beta_est[order(as.numeric(gsub(Mprefix, "", names(final_beta_est))))]
  final_beta_var <- unlist(final_beta_var)
  final_beta_var <- final_beta_var[order(as.numeric(gsub(Mprefix, "", names(final_beta_var))))]
  final_beta_pval <- pnorm(final_beta_est/sqrt(final_beta_var), lower.tail = F)
  final_beta_fdr = p.adjust(final_beta_pval, method = "fdr")
  final_sig_beta = final_beta_fdr < 0.05
  # map taxa names back to their original names
  name_map_vector <- setNames(MicrobNameMap[1,], MicrobNameMap[2,])
  names(final_beta_est) <- name_map_vector[names(final_beta_est)]
  names(final_beta_var) <- name_map_vector[names(final_beta_var)]
  cat("Final Averaging ended at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  flush.console()  # Ensures it prints immediately in interactive environments
  
  
  return(list(results = results,
              results_ref = results_ref,
              final_results = final_results,
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

# parallel computing inner function
parallel_job_cv <- function(
    g, # which group
    grouptaxaInd.lst, # stores taxa indices for all groups
    RefInd, # which reference in this group
    grouprefInd.lst, # stores reference taxa for all groups
    etaInd, # which eta candidate
    etaCand, # all eta candidates (vector)
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
    alpha_phi_results <- alpha_phi_Update_iloop_NR(gMicrobData = gMicrobData_train, CovDataWithIntcp = CovDataWithIntcp_train, eta = eta,
                                                   RefTaxonName = RefTaxonName, printloss = FALSE)
    # alpha.hat will be calculated using OLS approach and phiMat.hat from N-R
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
              group = g,
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
    alpha_phi_results <- alpha_phi_Update_iloop_NR(gMicrobData = gMicrobData_train, CovDataWithIntcp = CovDataWithIntcp_train, eta = eta,
                                                   RefTaxonName = RefTaxonName, printloss = FALSE)
    # alpha.hat will be calculated using OLS approach and phiMat.hat from N-R
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
  pairs_ave_SE = pairs_SE_sum/(pairs_count + (nrow(org.gMicrobData_val) - paris_count)*mech_prob)/phiMat.hat

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
    Newton_Raphson=FALSE,
    QRdecomp=FALSE
){
  
  # parallel computing: small jobs: (taxa group g, ref taxon r_g, eta candidate)
  gMicrobData <- MicrobData[,grouptaxaInd.lst[[g]]]
  
  eta <- etaCand[etaInd]
  
  # prepare dataset for later analysis
  CovDataWithIntcp <- cbind(1, CovData)
  testCovInd <- 1 + which(colnames(CovData)==testCovNewName)
  RefTaxonName <- colnames(MicrobData)[grouprefInd.lst[[g]]][RefInd]
  
  if (Newton_Raphson) {
    alpha_phi_results <- alpha_phi_Update_iloop_NR(gMicrobData = gMicrobData, CovDataWithIntcp = CovDataWithIntcp, eta = eta,
                                                   RefTaxonName = RefTaxonName, printloss = FALSE)
    # alpha.hat will be calculated using OLS approach and phiMat.hat from N-R
  } else {
    # OLS update
    alpha_phi_results <- alpha_phi_Update_iloop_OLS(gMicrobData = gMicrobData, CovDataWithIntcp = CovDataWithIntcp, eta = eta,
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
  
  
  return(list(eta = eta,
              group = g,
              # selection = selection,
              group_taxa_names = colnames(gMicrobData),
              ref_taxon_name = RefTaxonName,
              alpha_OLS = alpha.hat,
              phi_hat = phiMat.hat,
              # alpha_JK = alpha.JK,
              # var_alpha_JK = var.JK,
              beta_OLS = beta.mat
              # beta_JK = beta.JK.mat,
              # beta_var = var.beta.mat,
              # beta_p_val = beta.p_val
  ))
  
  
}

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
    jackknife_block = 10
){
  
  # parallel computing: small jobs: (taxa group g, ref taxon r_g, eta candidate)
  gMicrobData <- MicrobData[,grouptaxaInd.lst[[g]]]
  
  eta <- etaCand[etaInd]
  
  # prepare dataset for later analysis
  CovDataWithIntcp <- cbind(1, CovData)
  testCovInd <- 1 + which(colnames(CovData)==testCovNewName)
  RefTaxonName <- colnames(MicrobData)[grouprefInd.lst[[g]]][RefInd]
  
  if (Newton_Raphson) {
    alpha_phi_results <- alpha_phi_Update_iloop_NR(gMicrobData = gMicrobData, CovDataWithIntcp = CovDataWithIntcp, eta = eta,
                                                   RefTaxonName = RefTaxonName, printloss = FALSE)
    # alpha.hat will be calculated using OLS approach and phiMat.hat from N-R
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
  alpha.JK.diff = sweep(alpha.JK.mat, 2, alpha.JK, FUN = "-")
  # alpha.JK.diff <- alpha.JK.mat - matrix(rep(alpha.JK, jackknife_block), nrow = jackknife_block, byrow = TRUE)
  var.JK <- (jackknife_block-1)/jackknife_block * t(alpha.JK.diff)%*%alpha.JK.diff
  
  ## select taxa dependent on covariates of interest: Zg
  var.beta.mat <- matrix(diag(var.JK), nrow = ncol(CovDataWithIntcp))[testCovInd,]
  var.beta.mat[var.beta.mat==0] <- 1e-20 # since 0 was observed in variance and caused NaN issue for setting1 scenario1 with seed77
  
  ## z-test: p-value
  beta.p_val <- pnorm(abs(beta.JK.mat/sqrt(var.beta.mat)), lower.tail = FALSE)
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
# data-preprocessing
metaData=function(MicrobData,CovData,linkIDname,testCov=NULL,ctrlCov=NULL, # testCov are covariates of interest, ctrlCov are confounders
                  testMany=T,ctrlMany=F,MZILN=F){
  results=list()
  
  if(length(linkIDname)==0){
    stop("linkIDname is missing.")
  }
  
  if(length(testCov)>0 | length(ctrlCov)>0){
    if(sum(c(testCov,ctrlCov)%in%colnames(CovData))!=length(c(testCov,ctrlCov))){
      stop("Error: some covariates are not available in the data.")
    }
  }
  
  if(sum(testCov%in%ctrlCov)>0){
    cat("Warnings: Variables appeared in both testCov list and ctrlCov list will be treated as testCov.","\n")
  }
  
  # read microbiome data
  if(is.matrix(MicrobData))MdataWithId=data.matrix(MicrobData)
  if(is.data.frame(MicrobData))MdataWithId=data.matrix(MicrobData)
  if(is.character(MicrobData)){
    nCharac=nchar(MicrobData)
    if(substr(MicrobData,(nCharac-2),nCharac)=="csv"){
      MdataWithId=data.matrix(read.csv(file=MicrobData,header=T,na.strings=c("","NA")))
    }
    if(substr(MicrobData,(nCharac-2),nCharac)=="tsv"){
      MdataWithId=data.matrix(read.table(file=MicrobData, sep='\t',header=T,na.strings=c("","NA")))
    }
  }
  
  if(length(colnames(MdataWithId))!=ncol(MdataWithId))
    stop("Microbiome data lack variable names.")
  
  missPropMData=sum(is.na(MdataWithId[,linkIDname]))/nrow(MdataWithId)
  if(missPropMData>0.8){
    cat("Warning: There are over 80% missing values for the linkId variable in the Microbiome data file. 
               Double check the data format.","\n")
  }
  
  # read covariate data
  if(is.matrix(CovData))CovarWithId=data.matrix(CovData)
  if(is.data.frame(CovData))CovarWithId=data.matrix(CovData)
  if(is.character(CovData)){
    nCharac=nchar(CovData)
    if(substr(CovData,(nCharac-2),nCharac)=="csv"){
      CovarWithId=data.matrix(read.csv(file=CovData,header=T,na.strings=c("","NA")))
    }
    if(substr(CovData,(nCharac-2),nCharac)=="tsv"){
      CovarWithId=data.matrix(read.table(file=CovData, sep='\t',header=T,na.strings=c("","NA")))
    }
  }
  
  if(length(colnames(CovarWithId))!=ncol(CovarWithId))
    stop("Covariate data lack variable names.")
  
  missPropCovData=sum(is.na(CovarWithId[,linkIDname]))/nrow(CovarWithId)
  if(missPropCovData>0.8){
    cat("Warning: There are over 80% missing values for the linkId variable in the covariates data file. 
               Double check the data format.","\n")
  }
  
  Covariates1=CovarWithId[,!colnames(CovarWithId)%in%linkIDname,drop=F]
  
  # determine testCov and ctrlCov
  if(length(testCov)==0){
    if(!testMany){
      stop("No covariates are specified for estimating associations of interest.")
    }else{
      cat("Associations are being estimated for all covariates since no covariates are specified for testCov.","\n")
      testCov=colnames(Covariates1)
    }
  }
  results$testCov=testCov
  
  ctrlCov=ctrlCov[!ctrlCov%in%testCov]
  
  xNames=colnames(Covariates1)
  rm(Covariates1)
  
  if(length(ctrlCov)==0 & ctrlMany){
    cat("No control covariates are specified, 
            all variables except testCov are considered as control covariates.","\n")
    ctrlCov=xNames[!xNames%in%testCov]
  }
  results$ctrlCov=ctrlCov[!ctrlCov%in%testCov]
  
  # merge data to remove missing
  CovarWithId1=CovarWithId[,c(linkIDname,testCov,ctrlCov)]
  
  allRawData=data.matrix(na.omit(merge(CovarWithId1,MdataWithId,by=linkIDname,all.x=F,all.y=F)))
  
  CovarWithId=allRawData[,(colnames(allRawData)%in%colnames(CovarWithId1)),drop=F]
  Covariates=CovarWithId[,!colnames(CovarWithId)%in%linkIDname,drop=F]
  rm(CovarWithId1)
  
  if(!is.numeric(Covariates[,testCov,drop=F])){
    stop("There are non-numeric variables in the covariates for association test.")
  }
  
  MdataWithId=allRawData[,(colnames(allRawData)%in%colnames(MdataWithId))]
  Mdata_raw=MdataWithId[,!colnames(MdataWithId)%in%linkIDname,drop=F]
  rm(allRawData)
  
  # check zero taxa and subjects with zero taxa reads
  numTaxaNoReads=sum(colSums(Mdata_raw)==0)
  if(numTaxaNoReads>0){
    Mdata_raw=Mdata_raw[,!(colSums(Mdata_raw)==0)]
    cat("There are",numTaxaNoReads,"taxa without any sequencing reads and 
          excluded from the analysis","\n")
  }
  rm(numTaxaNoReads)
  
  numSubNoReads=sum(rowSums(Mdata_raw)==0)
  if(numSubNoReads>0){
    cat("There are",numSubNoReads,"subjects without any sequencing reads and 
          excluded from the analysis","\n")
    subKeep=!(rowSums(Mdata_raw)==0)
    Mdata_raw=Mdata_raw[subKeep,]
    MdataWithId=MdataWithId[subKeep,]
    rm(subKeep)
  }
  rm(numSubNoReads)
  
  Mdata=Mdata_raw
  
  rm(Mdata_raw)
  
  microbName=colnames(Mdata)
  newMicrobNames=paste0("microb",seq(length(microbName)))
  results$Mprefix="microb"
  
  colnames(Mdata)=newMicrobNames
  
  MdataWithId_new=cbind(MdataWithId[,linkIDname,drop=F],Mdata)
  results$microbName=microbName
  results$newMicrobNames=newMicrobNames
  rm(microbName,newMicrobNames)
  
  xNames=colnames(Covariates)
  nCov=length(xNames)
  
  if(sum(is.na(Covariates))>0){
    cat("Samples with missing covariate values are removed from the analysis.","\n")
  }
  
  if(!is.numeric(Covariates[,ctrlCov,drop=F])){
    cat("Warnings: there are non-numeric variables in the control covariates","\n")
    nCtrlCov=length(ctrlCov)
    numCheck=unlist(lapply(seq(nCtrlCov),function(i)is.numeric(Covariates[,ctrlCov[i]])))+0
    for(i in which(numCheck==0)){
      Covariates[,ctrlCov[i]]=as.numeric(factor(Covariates[,ctrlCov[i]]))
    }
  }
  
  colnames(Covariates)=xNames
  binCheck=unlist(lapply(seq(nCov),function(i)dim(table(Covariates[,xNames[i]]))))
  
  if(length(which(binCheck==2))>0){
    Covariates=Covariates[,c(xNames[binCheck!=2],xNames[binCheck==2]),drop=F]
    binaryInd=length(which(binCheck!=2))+1
    results$varNamForBin=xNames[binCheck==2]
    results$BinVars=length(results$varNamForBin)
    for(i in results$varNamForBin){
      mini=min(Covariates[,i],na.rm=T)
      maxi=max(Covariates[,i],na.rm=T)
      if(!(mini==0 & maxi==1)){
        Covariates[Covariates[,i]==mini,i]=0
        Covariates[Covariates[,i]==maxi,i]=1
        cat("Binary covariate",i,"is not coded as 0/1 which may generate analysis bias. It has been changed to 0/1. The changed covariates data can be extracted from the result file.","\n")
      }
    }
    #cat(length(which(binCheck==2)),"binary covariates are detected.","\n")
  }else{
    results$BinVars=0
    binaryInd=NULL
    results$varNamForBin=NULL
  }
  
  # find the position of the first binary predictor and the rest are all binary
  results$binaryInd=binaryInd
  results$xNames=colnames(Covariates)  
  xNewNames=paste0("x",seq(length(xNames)))
  colnames(Covariates)=xNewNames
  results$covsPrefix="x"
  results$xNewNames=xNewNames
  
  results$testCovInd=which((results$xNames)%in%testCov)
  results$testCovInOrder=results$xNames[results$testCovInd]
  results$testCovInNewNam=results$xNewNames[results$testCovInd]
  rm(xNames,xNewNames)
  
  CovarWithId_new=cbind(CovarWithId[,linkIDname,drop=F],Covariates)
  
  data=merge(MdataWithId_new, CovarWithId_new,by=linkIDname,all.x=F,all.y=F)
  results$covariatesData=CovarWithId_new
  colnames(results$covariatesData)=c(linkIDname,results$xNames)
  rm(MdataWithId_new,CovarWithId_new)
  results$data=na.omit(data)
  rm(data)
  cat("Data dimensions (after removing missing data if any):","\n")
  cat(dim(results$data)[1],"samples","\n")
  cat(ncol(Mdata),"OTU's or microbial taxa","\n")
  
  if(!MZILN)cat(length(results$testCovInOrder),"testCov variables in the analysis","\n")
  if(MZILN)cat(length(results$testCovInOrder),"covariates in the analysis","\n")
  
  if(length(results$testCovInOrder)>0){
    if(!MZILN)print("These are the testCov variables:")
    if(MZILN)print("These are the covariates:")
    print(testCov)
  }
  rm(testCov)
  if(!MZILN){
    cat(length(results$ctrlCov),"ctrlCov variables in the analysis ","\n")
    if(length(results$ctrlCov)>0){
      print("These are the ctrlCov variables:")
      print(ctrlCov)
    }
    rm(ctrlCov)
  }
  cat(results$BinVars,"binary covariates in the analysis","\n")
  if(results$BinVars>0){
    print("These are the binary covariates:")
    print(results$varNamForBin)
  }
  rm(Mdata,Covariates,binCheck)
  return(results)
}



# calculate eta given p
estiEta=function(
    MicrobData,
    mechZeroPerc,
    a=10^(-2),
    b=10
){
  numeData=c(as.numeric(data.matrix(MicrobData)))
  nDataPoints=length(numeData)
  allValues=unique(sort(numeData)) # some extreme values
  
  nDisticValues=length(allValues)
  allZeroPerc=sum(numeData==0)/nDataPoints
  
  postiveValues=allValues[-1]
  
  empProps=rep(NA,(nDisticValues-1))
  for(j in 1:(nDisticValues-1)){
    empProps[j]=sum(numeData==(postiveValues[j]))/nDataPoints
  }
  
  fEta=function(x){
    equat=sum(empProps*exp(-x*postiveValues)/(1-exp(-x*postiveValues))) - mechZeroPerc*allZeroPerc
    return(equat)
  }
  
  solveEqu=uniroot(f=fEta,interval=c(a,b),extendInt="downX")
  
  eta=solveEqu$root
  rm(solveEqu)
  return(eta)
}

# replace Inf by a finite large value
InfCorrect <- function(para, correct_scale = 1e+30) {
  para[which(is.infinite(para))] <- correct_scale*sign(para[which(is.infinite(para))])
  return(para)
}

# QR decomposition: Ax = y, A is ZZsum, x is alpha to solve, y is ZUsum
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

# OLS to update alpha and phi iteratively
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
      calcloss_iloop(phiMat, k2k1.lst, W.lst, U.lst, AZ.lst, twoPosSub, alpha.new)
    }
  }
  
  
  # store all updates of alpha and phi for later debugging / check convergence performance
  results$alpha_update_mat <- alpha_update.mat
  results$phi_update_mat <- phi_update.mat
  
  
  return(results)
}


# # Newton-Raphson to update alpha and phi iteratively
# ## update alpha^{K+1} and phi: outer loop k1k2
# alpha_phi_Update_iloop_NR <- function(
    #     gMicrobData, # N by (K+1) taxa matrix
#     CovDataWithIntcp, # N by (1+Q+S) covariate matrix which includes intercept
#     eta, # parameter of zero probability mechanism
#     RefTaxonName, 
#     alpha.init = rep(0, 4*(ncol(gMicrobData)-1)), # initial guess of coefficients
#     phiMat.init = matrix(1, nrow = ncol(gMicrobData), ncol = ncol(gMicrobData)),
#     printloss = FALSE,
#     abnormal_ginv_threshold = 0.025
# ){
#   results <- list()
#   force(alpha.init)  # Ensures evaluation before parallel execution
#   force(phiMat.init)  # Ensures evaluation before parallel execution
#   
#   # calculate K: K+1 is the number of taxa
#   K <- ncol(gMicrobData) - 1
#   
#   # reorganize the microbiome matrix by moving the reference taxon to the last column
#   refTaxonInd <- which(colnames(gMicrobData)==RefTaxonName)
#   org.gMicrobData <- cbind(gMicrobData[,-refTaxonInd], gMicrobData[,refTaxonInd])
#   rm(gMicrobData)
#   
#   # to calculate and store AZ and Weight matrices which are not dependent on alpha --> CANNOT avoid calculating AZWAZ or AZW them in every while loop
#   k2k1.lst <- list()
#   AZ.lst <- list()
#   W.lst <- list()
#   U.lst <- list()
#   twoPos.vec <- rep(FALSE, nrow(org.gMicrobData))
#   for (i in 1:nrow(org.gMicrobData)) {
#     # extract positions of positive taxa
#     taxa.nonzero.pos=which(org.gMicrobData[i,]!=0)
#     # only consider subjects with >= 2 positive taxa
#     if (length(taxa.nonzero.pos) >= 2) {
#       twoPos.vec[i] <- TRUE
#       pairs=combn(x=sort(taxa.nonzero.pos,decreasing=T),m=2)
#       k2k1.lst[[i]] <- pairs
#       
#       # create Ai, then calculate AZi
#       Ai <- matrix(0, nrow = ncol(pairs), ncol = K)
#       Ai[cbind(1:ncol(pairs), pairs[2,])] <- 1 # (l,kl)
#       Ai[cbind((1:ncol(pairs))[pairs[1,] < (K+1)], pairs[1,(pairs[1,] < (K+1))])] <- -1
#       AZ.lst[[i]] <- t(kronecker(t(Ai), t(CovDataWithIntcp[i,])))
#       
#       # Wi
#       W.lst[[i]] <- as.numeric(1/sqrt((1-exp(-eta*(org.gMicrobData[i,pairs[2,]])))*(1-exp(-eta*(org.gMicrobData[i,pairs[1,]])))))
#       
#       # Ui
#       U.lst[[i]] <- as.numeric(log(org.gMicrobData[i,pairs[2,]] / org.gMicrobData[i,pairs[1,]]))
#     }
#   }
#   twoPosSub <- which(twoPos.vec)
#   
#   # store parameters that won't change in loop
#   results$k2k1.lst <- k2k1.lst
#   results$AZ.lst <- AZ.lst
#   results$W.lst <- W.lst
#   results$U.lst <- U.lst
#   results$twoPosSub <- twoPosSub
#   
#   
#   # update alpha and phi
#   ## alpha is the coefficient
#   alpha_update.mat <- matrix(alpha.init, nrow = 1) # to store all updates of alpha
#   phi_update.mat <- matrix(as.vector(t(phiMat.init)), nrow = 1)
#   
#   alpha.new <- alpha.init
#   alpha.old <- rep(-100, 4*K)
#   phiMat = phiMat.init
#   
#   # store parameters which are updated in loop
#   results$phiMat <- phiMat
#   results$alpha.hat <- as.numeric(alpha.new) # depends on phi
#   
#   
#   # start the updating loop
#   ite = 0
#   while (mean(abs(alpha.new - alpha.old)) > 1e-3 & ite < 20) {
#     ite <- ite + 1
#     alpha.old <- alpha.new
#     
#     # initialize Score, Hessian, and phi's numerator&denominator to zeros for this new loop of update
#     S <- rep(0, 4*K)
#     H <- matrix(0, nrow = 4*K, ncol = 4*K)
#     phi_nume <- matrix(0, nrow = K+1, ncol = K+1)
#     phi_deno <- matrix(0.01, nrow = K+1, ncol = K+1)
#     phi_deno[upper.tri(phi_deno)] <- 0
#     
#     # only use subjects who have >= 2 positive taxa
#     for (i in twoPosSub) {
#       pairs <- k2k1.lst[[i]]
#       phi_i <- phiMat[cbind(pairs[2,], pairs[1,])]
#       AZi <- AZ.lst[[i]]
#       AZalpha <- AZi%*%alpha.old
#       Wi <- W.lst[[i]]
#       Ui <- U.lst[[i]]
#       
#       # update Hessian and Score
#       WoverPhi <- InfCorrect(Wi/phi_i)
#       AZW <- InfCorrect(t(AZi)%*%diag(WoverPhi))
#       
#       Hi <- InfCorrect(AZW%*%AZi)
#       H <- InfCorrect(H + Hi)
#       
#       U_AZalpha <- InfCorrect(Ui - AZalpha)
#       Si <- InfCorrect(AZW%*%U_AZalpha)
#       S <- InfCorrect(S + Si)
#       
#       # update phi's numerator and denominator
#       weighted.squares <- InfCorrect(Wi*(Ui - AZalpha)^2)
#       phi_nume[cbind(pairs[2,], pairs[1,])] <- InfCorrect(phi_nume[cbind(pairs[2,], pairs[1,])] + weighted.squares)
#       phi_deno[cbind(pairs[2,], pairs[1,])] <- InfCorrect(phi_deno[cbind(pairs[2,], pairs[1,])] + Wi)
#       
#     }
#     # update phi
#     phiMat <- phi_nume / phi_deno
#     
#     
#     # if the pseudo-inverse cannot be properly calculated: break this while-loop and use the previous alpha.hat as final result
#     H.inv <- ginv(as.matrix(H), tol = 1e-40)
#     sum_diag_identity <- sum(abs(diag(H.inv%*%H)))
#     if (abs(sum_diag_identity - nrow(H))/nrow(H) > abnormal_ginv_threshold) break
#     # else: there is nothing wrong with ginv(ZZsum), update alpha
#     alpha.new <- alpha.old + H.inv%*%S
#     phi_update.mat <- rbind(phi_update.mat, as.vector(t(phiMat)))
#     alpha_update.mat <- rbind(alpha_update.mat, as.numeric(alpha.new))
#     
#     
#     # store the latest parameter values
#     results$phiMat <- phiMat
#     results$alpha.hat <- as.numeric(alpha.new) # depends on phi
#     
#     if (printloss) {
#       loss <- calcloss_iloop(phiMat, k2k1.lst, W.lst, U.lst, AZ.lst, twoPosSub, alpha.new)
#     }
#   }
#   
#   
#   # store all updates of alpha and phi for later debugging / check convergence performance
#   results$alpha_update_mat <- alpha_update.mat
#   results$phi_update_mat <- phi_update.mat
#   
#   
#   return(results)
# }



# ## loss function
# calcloss_iloop <- function(
    #     phiMat,
#     k2k1.lst,
#     W.lst,
#     U.lst,
#     AZ.lst,
#     twoPosSub, # vector of indices of which subjects have >= 2 positive taxa
#     alpha.hat
# ){
#   loss <- 0
#   for (i in twoPosSub) {
#     pairs <- k2k1.lst[[i]]
#     phi_i <- phiMat[cbind(pairs[2,], pairs[1,])]
#     AZi <- AZ.lst[[i]]
#     AZalpha <- AZi%*%alpha.hat
#     Wi <- W.lst[[i]]
#     Ui <- U.lst[[i]]
#     
#     loss <- loss + sum((Wi/phi_i)*(Ui - AZalpha)^2)
#   }
#   
#   return(loss)
# }



# # get true coefficients
# getTrueAlpha <- function(
    #     gMicrobData_microbnames,
#     RefTaxonName,
#     coefMat
# ){
#   nTaxa <- length(gMicrobData_microbnames)
#   microbInd <- as.numeric(gsub("microb", "", gMicrobData_microbnames))
#   beta.mat <- coefMat[,microbInd]
#   colnames(beta.mat) <- gMicrobData_microbnames
#   
#   refTaxonInd <- which(gMicrobData_microbnames==RefTaxonName)
#   org.beta.mat <- cbind(beta.mat[,-refTaxonInd], beta.mat[,refTaxonInd])
#   rm(beta.mat)
#   
#   trueAlpha <- org.beta.mat[,1:(nTaxa-1)] - matrix(rep(org.beta.mat[,nTaxa], nTaxa-1), ncol = nTaxa-1)
#   
#   return(trueAlpha)
# }



# # get true phi
# getTruePhi <- function(
    #     gMicrobData_microbnames,
#     RefTaxonName,
#     covNorm
# ){
#   nTaxa <- length(gMicrobData_microbnames)
#   refTaxonInd <- which(gMicrobData_microbnames==RefTaxonName)
#   org.gMicrobData_microbnames <- c(gMicrobData_microbnames[-refTaxonInd], gMicrobData_microbnames[refTaxonInd])
#   microbInd <- as.numeric(gsub("microb", "", org.gMicrobData_microbnames))
#   
#   # extract the submatrix
#   gCovMat <- covNorm[microbInd, microbInd]
#   rm(covNorm)
#   
#   # Phi matrix: upper-right triangle matrix
#   truePhiMat <- matrix(0, nrow = nTaxa, ncol = nTaxa)
#   for (k1 in 1:(nTaxa-1)) {
#     for (k2 in (k1+1):nTaxa) {
#       truePhiMat[k1,k2] <- gCovMat[k1,k1] + gCovMat[k2,k2] - 2*gCovMat[k1,k2]
#     }
#   }
#   
#   return()
# }

