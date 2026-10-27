#' IWEE: Inverse-Weighted Estimating Equations for Zero-Inflated Microbiome Data
#'
#' Implements an inverse-weighted estimating equation (IWEE) approach for
#' identifying taxa associated with covariates of interest in zero-inflated,
#' compositional microbiome data. The main function is [IWEE()]; [metaData()]
#' and [estiEta()] expose the data preparation and the zero-mechanism
#' parameter calculation, and [load_iwee_example()] provides example data.
#'
#' @keywords internal
#' @import Matrix
#' @importFrom MASS ginv
#' @importFrom foreach foreach %dopar% %:% getDoParWorkers
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom glmnet glmnet
#' @importFrom stats p.adjust pnorm sd setNames
#' @importFrom utils combn flush.console
"_PACKAGE"

# Loop variables created by foreach() and picked up by the %dopar% bodies
utils::globalVariables(c("RefInd", "etaInd", "val_fold", "g"))
