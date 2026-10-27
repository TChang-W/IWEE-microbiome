#' Load the example dataset shipped with the package
#'
#' Returns a simulated microbiome dataset used in the examples and tests. The
#' data were generated with 100 samples and 400 taxa, a binary covariate of
#' interest `X1` and two confounders `W1` and `W2`; about 42% of the
#' abundance entries are zero. The file is stored in `inst/extdata` and was
#' produced by the simulation file `paras_nTaxa400_1_binary_X_0.94_24.12_26.rds`
#' (number of taxa 400, random seed 1, binary `X`, signal-to-noise ratio 0.94,
#' floor-zero percentage 24.12).
#'
#' @return A named list with the elements
#' \describe{
#'   \item{`MicrobData`}{Numeric matrix (100 by 401) of observed abundances
#'     with the identifier column `id` followed by `taxon1` to `taxon400`.}
#'   \item{`CovData`}{Numeric matrix (100 by 4) with columns `id`, `X1`
#'     (covariate of interest), `W1` and `W2` (confounders).}
#'   \item{`trueZeroAA`, `sampAA`}{Matrices (100 by 401, with `id` column) of
#'     intermediate simulation outputs with the same layout as `MicrobData`.}
#'   \item{`coefMat`, `truBeta`}{Coefficient matrices (4 by 400) with rows
#'     `intercept`, `X1`, `W1`, `W2` and columns `taxon1` to `taxon400`; the
#'     `X1` row of `truBeta` holds the true associations to recover.}
#'   \item{`covNorm`}{400 by 400 matrix used in the simulation.}
#'   \item{`true_eta`}{True value of \eqn{\eta} (0.45).}
#'   \item{`Ci_rate`, `Ci_setting`, `trueZeroPerc`, `falseZeroPerc`,
#'     `floorZeroPerc`, `mechanismZeroPerc`, `totalZeroPerc`, `SNR`}{Scalars
#'     describing the simulation settings and the resulting zero
#'     percentages.}
#' }
#'
#' @examples
#' dat <- load_iwee_example()
#' dim(dat$MicrobData)
#' head(dat$CovData)
#' # True non-zero effects of X1
#' which(dat$truBeta["X1", ] != 0)
#' @export
load_iwee_example <- function() {
  path <- system.file("extdata", "iwee_example.rds", package = "IWEE")
  if (!nzchar(path)) stop("The example data file was not found in the IWEE installation.")
  readRDS(path)
}
