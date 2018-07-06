#' speedranger
#'
#' @name speedranger
#' @docType package
#' @useDynLib speedranger
#' @importFrom Rcpp sourceCpp
#' @importFrom parallel mclapply detectCores
#' @importFrom Rsamtools BamFile
#' @importFrom data.table data.table setkey as.data.table fread fwrite setDT rbindlist dcast.data.table copy
#' @importFrom stringr str_sub str_match
NULL
