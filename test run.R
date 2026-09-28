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
source("IWEE_utils_withoutSparsematrix_QRdecomp_v3.R")
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
source("IWEE_utils_withoutSparsematrix_QRdecomp_v3.R")
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

