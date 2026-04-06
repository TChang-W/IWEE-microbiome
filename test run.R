library(MASS) # ginv()
library(Matrix) # for QR decomposition
library(foreach)
library(doParallel)
library(peakRAM)
library(lmerTest)

# paras_213_binary_X <- readRDS("paras/nTaxa50/paras_nTaxa50_213_binary_X_3.47.rds")
# paras_213_binary_X <- readRDS("paras/nTaxa2/paras_nTaxa2_213_binary_X_1.16.rds")
paras_213_binary_X <- readRDS("con_X/paras/diffzero/nTaxa400_floor30/paras_nTaxa400_213_con_X_10.00_14.79_8.rds")

paras_213_binary_X <- readRDS("bin_X-Copy/paras/diffzero/nTaxa400_floor30/paras_nTaxa400_213_binary_X_0.91_11.61_3.rds")
paras_213_binary_X = paras_nTaxa400_213_binary_X_0.91_11.61_3

list2env(paras_213_binary_X, envir = .GlobalEnv)
source("IWEE_utils_withoutSparsematrix_QRdecomp_v2.R")
iwee_res <- IWEE(
    MicrobData,
    CovData,
    num_cores = 12, 
    testCov = "X1", # covariates of interest
    # ctrlCov = paste0("W", 1:2), # confounders
    ctrlCov = NULL, # confounders
    linkIDname = "id",
    QRdecomp = TRUE,
    nGroup = 2,
    nEta = 11, # number of eta candidates corresponding to p = 0 to 1
    Ng = 24, # the smallest number of taxa in a taxa group after partitioning
    nRef = 20, # the number of selected reference candidates in each group
    nBestRef = 2, # the number of final selected best reference taxa in each group
    k_fold = 3, # number of folds 
    testMany = F,
    ctrlMany = T,
    randseed = 213
)
filename = sprintf("./sim1_ci2_results/iweeres_%d_%s_%.2f.rds", randseed, "binary_X", SNR)
saveRDS(iwee_res, file = filename)


num_cores = 12
testCov = "X1"
ctrlCov = paste0("W", 1:2)
linkIDname = "id"
QRdecomp = TRUE
nGroup = 3
nEta = 11 # number of eta candidates corresponding to p = 0 to 1
Ng = 24 # the smallest number of taxa in a taxa group after partitioning
nRef = 20 # the number of selected reference candidates in each group
nBestRef = 2 # the number of final selected best reference taxa in each group
k_fold = 3 # number of folds 
testMany = F
ctrlMany = T
randseed = 213


paras_213_binary_X <- readRDS("paras/nTaxa3/paras_nTaxa3_213_binary_X_0.76.rds")

paras_213_binary_X <- readRDS("paras/nTaxa2/paras_nTaxa2_213_binary_X_0.84.rds")
list2env(paras_213_binary_X, envir = .GlobalEnv)
source("IWEE_utils_withoutSparsematrix_QRdecomp_v2.R")
iwee_res <- IWEE(
  MicrobData,
  CovData,
  num_cores = 12, 
  testCov = "X1", # covariates of interest
  ctrlCov = paste0("W", 1:2), # confounders
  linkIDname = "id",
  QRdecomp = TRUE,
  nGroup = 1,
  nEta = 11, # number of eta candidates corresponding to p = 0 to 1
  Ng = 2, # the smallest number of taxa in a taxa group after partitioning
  nRef = 2, # the number of selected reference candidates in each group
  nBestRef = 1, # the number of final selected best reference taxa in each group
  k_fold = 3, # number of folds 
  testMany = F,
  ctrlMany = T,
  randseed = 213
)
filename = sprintf("./sim1_ci2_results/iweeres_%d_%s_%.2f.rds", randseed, "binary_X", SNR)
saveRDS(iwee_res, file = filename)



############### lm test ###############

MicrobData_tmp = data.frame(cbind(MicrobData, CovData))
MicrobData_tmp$U = log(MicrobData_tmp$microb2 / MicrobData_tmp$microb1) 
MicrobData_tmp2 = subset(MicrobData_tmp, is.finite(U))
lm_model = lm(U ~ x1, data = MicrobData_tmp2)
lm_model

MicrobData_tmp2$W = sqrt(1/((1-exp(-eta*MicrobData_tmp2$microb1))*(1-exp(-eta*MicrobData_tmp2$microb2))))
MicrobData_tmp2  = MicrobData_tmp2 %>%
  mutate(XW = x1*W,
         UW = U*W)
lm_model2 = lm(UW ~ -1 + XW + W, data = MicrobData_tmp2)
lm_model2

sampAA = data.frame(sampAA)
sampAA$U = log(sampAA$taxon2 / sampAA$taxon1) 
MicrobData = trueZeroAA
MicrobData_tmp = data.frame(cbind(MicrobData, CovData))
MicrobData_tmp2 = subset(MicrobData_tmp, is.finite(sampAA$U))
MicrobData_tmp2$U = log(MicrobData_tmp2$taxon2 / MicrobData_tmp2$taxon1) 
lm_model = lm(U ~ X1, data = MicrobData_tmp2)



############### SE test #################
val_fold = 1
g = 1
etaInd = 5
eta = etaCand[etaInd]
RefInd = 1
RefTaxonName <- colnames(MicrobData)[grouprefInd.lst[[g]]][RefInd]

gMicrobData <- MicrobData[,grouptaxaInd.lst[[g]]]
test_index = groupfoldInd.lst[[val_fold]]

# SE from the hand 
fold_alpha = results[[1]][[1]][[etaInd]][[val_fold]]$alpha_OLS
fold_phi = results[[1]][[1]][[etaInd]][[val_fold]]$alpha_phi_results$phiMat[1, 1]

MicrobData_tmp = data.frame(cbind(MicrobData, CovData))
MicrobData_tmp$U = log(MicrobData_tmp$microb2 / MicrobData_tmp$microb1) 
microb_error = (MicrobData_tmp$U - MicrobData_tmp$x1*fold_alpha[2] - fold_alpha[1])^2/fold_phi

microb_error_val = microb_error[test_index]
microb_error_val = microb_error_val[is.finite(microb_error_val)]

# validation
gMicrobData_val = gMicrobData[test_index, ]
CovDataWithIntcp_val = CovDataWithIntcp[test_index, ]
# training
gMicrobData_train = gMicrobData[-test_index, ]
CovDataWithIntcp_train = CovDataWithIntcp[-test_index, ]

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

org.gMicrobData_val <- cbind(gMicrobData_val[,-RefInd], gMicrobData_val[,RefInd])
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
    
    ####!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    # AZi%*%alpha.hat should be a scalar, return matrix instead
    phi_i <- phiMat.hat[cbind(pairs[2,], pairs[1,])]
    SE_eta[i] <- mean((Ui - AZi%*%alpha.hat)^2/phi_i, na.rm = TRUE)
    
    # ignore phi_i
    # SE_eta[i] <- mean((Ui - AZi%*%alpha.hat)^2, na.rm = TRUE)
  }else{
    SE_eta[i] = NA
  }
}

mean(microb_error_val^2/fold1_phi, na.rm = T)
