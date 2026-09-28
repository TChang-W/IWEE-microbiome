library(MASS)

sigmoid = function(x) {
  1 / (1 + exp(-x))
}

# true Beta for X
set.seed(3)
nTaxa = 400
N = 100
truBeta = matrix(0, 4, nTaxa)
truBeta[1, ] = runif(nTaxa, min = 10, max = 10)

#################################### params ####################################

# Ci_mat = matrix(c(1/2,1/2, 1/2,1/6, 1/2,1/10, 1/2,1/14, 1/2,1/20), nrow=2)
# Ci_mat = matrix(c(1/300,1/300, 1/2,1/6, 1/2,1/10, 1/2,1/14, 1/2,1/20), nrow=2)
Ci_max_lst = 1/c(seq(100, 1000, by = 100), seq(2000, 10000, by = 1000))
Ci_min_lst = 1/c(seq(100, 1000, by = 100), seq(2000, 10000, by = 1000))
Ci_scale = 100

para = 60
b_setting = 10

i = 1
j = 19
eta = 0.1

Ci_min = Ci_min_lst[j]
Ci_max = Ci_max_lst[i]
u = runif(nTaxa)
# SNR = 1.05
truBeta[2, ] = ifelse(u < 0.75, 0, ifelse(u > 0.9, runif(nTaxa, 3, 4), runif(nTaxa, 4, 5)))
covNorm <- diag(runif(nTaxa, 2, 3))
randeff_var = 0.5
# covNorm <- diag(rep(3.5, nTaxa))
# randeff_var = runif(1, 0.25, 0.75)

# SNR = 3.47
# truBeta[2, ] = ifelse(u < 0.75, 0, ifelse(u > 0.9, runif(nTaxa, 2, 3), runif(nTaxa, 4, 5)))
# covNorm <- diag(runif(nTaxa, 0.1, 0.2))
# randeff_var = 0.5

# # SNR = 5.54
# truBeta[2, ] = ifelse(u < 0.75, 0, ifelse(u > 0.9, runif(nTaxa, 2, 3), runif(nTaxa, 4, 5)))
# covNorm <- diag(runif(nTaxa, 0.05, 0.1))
# randeff_var = 0.1

# # SNR = 8
# truBeta[2, ] = ifelse(u < 0.75, 0, ifelse(u > 0.9, runif(nTaxa, 3, 4), runif(nTaxa, 3, 4)))
# covNorm <- diag(runif(nTaxa, 4, 5))
# randeff_var = 0.5

################################################################################
# truBeta[3, ] = runif(nTaxa, min = -3, max = 3)
# truBeta[4, ] = runif(nTaxa, min = -3, max = 3)
truBeta[3, ] = 0
truBeta[4, ] = 0

rownames(truBeta) = c("intercept", "X1", "W1", "W2")
colnames(truBeta) = paste0("taxon", seq_len(nTaxa))
write.table(truBeta, file = "./simTrueBetaMat_con.csv", sep = ",")

if (nTaxa > ncol(truBeta)) stop(sprintf("The number of taxa cannot be greater than %d.", ncol(truBeta)))
coefMat <- truBeta[,1:(nTaxa)]

# generate covariates and confounders
n_omega <- nrow(truBeta)-1 # omega1, omega2, omega3
Sigma <- outer(1:n_omega, 1:n_omega, function(i, j) rho^abs(i - j))

# remove W1, W2
# truBeta = truBeta[1:2, ]
# coefMat = coefMat[1:2, ]


for (i in 1:50) {
  
  omega.mat <- mvrnorm(N, mu = rep(0, n_omega), Sigma = Sigma)
  CovData <- cbind(runif(N, min = -0.5, max = 1), as.numeric(omega.mat[,2]>0), omega.mat[,3])
  # CovData = CovData[, 1, drop = FALSE]
  CovDataWithIntcp <- cbind(1, CovData)
  
  CovData <- cbind(c(1:N), CovData)
  colnames(CovData) <- c("id", "X1", paste0("W", 1:2))
  # colnames(CovData) <- c("id", "X1")
  
  
  # generate the true \mathcal{Yi} (in gut) from log-normal distribution
  ## generate mean for mulrivariate normal distribution
  muNorm <- CovDataWithIntcp%*%coefMat
  randeff <- rnorm(N, sd = sqrt(randeff_var))
  ## generate multivariate normal residual error
  ### temporarily assume independent taxa
  
  epsMat <- mvrnorm(N, mu = rep(0, nTaxa), Sigma = covNorm)
  ## generate the true \mathcal{Yi} (in gut)
  ecoAA <- 2^(muNorm + randeff + epsMat)
  var_error <- diag(var(randeff + epsMat))
  
  
  # Sparsity
  ## true zero 
  z.mat <- mvrnorm(N, mu = rep(0, nTaxa), Sigma = diag(1, nrow=nTaxa, ncol=nTaxa))
  cp <- qnorm(true_zero_sparsity)
  trueZero.mat <- (z.mat > cp)
  trueZero.mat <- matrix(as.numeric(trueZero.mat), nrow = nrow(trueZero.mat))
  # (no true zero)
  trueZero.mat = matrix(1, nrow = N, ncol = nTaxa)
  ## observed AA with true zero
  trueZeroAA <- ecoAA*trueZero.mat
  trueZeroPerc <- 100*mean(trueZeroAA==0)
  
  
  # Mimic sequencing
  ## generate sequencing depth Ci depending on X1
  b <- b_setting
  Ci <- pmax(Ci_min, pmin(sigmoid(b_setting * CovData[,2])/Ci_scale, Ci_max))
  # plot(CovData[,2], Ci)
  ## generate sample truth Yi
  sampAA <- Ci*trueZeroAA
  sampAA <- floor(sampAA)
  floorZeroPerc <- 100*mean(sampAA==0) - trueZeroPerc
  
  ## true and false zero (probability mechanism)
  ZeroProb <- exp(-eta*sampAA) # this mechanism include both the true and false zeros
  Zero.mat <- 1 - matrix(rbinom(length(ZeroProb), size=1, prob=ZeroProb), nrow=nrow(ZeroProb))
  # Zero.mat = 1
  # Scenario 4: eta=0.07,FZPerc=49%; eta=0.2,FZP=30%; eta=0.55,FZP=20%
  # eta=0.07: Scenario4,FZPerc=49.08%; 
  mechanismZeroPerc <- 100*(mean(Zero.mat==0) - mean(sampAA==0))
  falseZeroPerc <- mechanismZeroPerc + floorZeroPerc
  
  # generate the observed AA
  MicrobData <- sampAA*Zero.mat
  # MicrobData <- sampAA
  # MicrobData <- ceiling(MicrobData) # avoid using floor which increases/disturb the sparsity
  
  # add ID column for both MicrobData and CovData
  MicrobData <- cbind(c(1:N), MicrobData)
  colnames(MicrobData) <- c("id", paste0("taxon", 1:(nTaxa)))
  
  
  ecoAA <- cbind(c(1:N), ecoAA)
  colnames(ecoAA) <- c("id", paste0("taxon", 1:(nTaxa)))
  
  
  sampAA <- cbind(c(1:N), sampAA)
  colnames(sampAA) <- c("id", paste0("taxon", 1:(nTaxa)))
  
  trueZeroAA <- cbind(c(1:N), trueZeroAA)
  colnames(trueZeroAA) <- c("id", paste0("taxon", 1:(nTaxa)))
  
  X_beta = coefMat[2, ]
  X_selection = X_beta!=0
  var_mu = X_beta * var(CovDataWithIntcp[, 2]) * X_beta
  
  SNR = mean(sqrt(var_mu/var_error)[X_selection]) # take beta != 0 subset, take average
  SNR = round(SNR, 2)

  # add SNR, beta, save everything into .RDS
  results <- list(trueZeroAA = trueZeroAA,
                  sampAA = sampAA,
                  MicrobData = MicrobData,
                  CovData = CovData,
                  coefMat = coefMat,
                  covNorm = covNorm,
                  b_setting = b_setting,
                  true_eta = eta,
                  truBeta = truBeta,
                  trueZeroPerc = trueZeroPerc,
                  falseZeroPerc = falseZeroPerc,
                  floorZeroPerc = floorZeroPerc,
                  mechanismZeroPerc = mechanismZeroPerc,
                  totalZeroPerc = trueZeroPerc + falseZeroPerc,
                  SNR = round(SNR, 2))

  tail(results, -6)
  
  # # floor zero filename
  # filename = sprintf("./con_X/paras/diffzero/nTaxa%d_floor%d/paras_nTaxa%d_%d_%s_%.2f_%.2f_%d.rds", nTaxa, para, nTaxa, randseed, "con_X", b_setting, para, i)
  
  # # sparsity filename
  # filename = sprintf("./con_X/paras/diffspars/nTaxa%d_spars%d/paras_nTaxa%d_%d_%s_%.2f_%.2f_%d.rds", nTaxa, para, nTaxa, randseed, "con_X", b_setting, para, i)
  
  
  saveRDS(results, file = filename)
  # return(results)
}
