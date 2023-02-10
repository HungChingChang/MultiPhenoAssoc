#' Subdataset from METABRIC
#' (https://www.mercuriolab.umassmed.edu/metabric)
#'
#' @description Containing 50 samples with grade 2 and 50 samples with grade 3.
#' 3k genes are pre-selected by coefficient of variation.
#'
#' @format A list with three types of data
#' \describe{
#'   \item{Pheno}{A matrix with four phenotypes.}
#'   \item{Gene.expr}{A matrix with 3k genes}
#'   \item{Confounder: }{A matrix with two columns. The first column is age, and the second column is generated from normal distribution with mean 20 and variance 10}
#' }
"exampleData"
