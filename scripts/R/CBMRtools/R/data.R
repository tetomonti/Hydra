#' Example dataset
#' TCGA breast cancer RNA-seq expression set (mad filtered 100 genes), log2 transformed
#' 
#' @docType data
#' @name eSet.brca.100
#' @usage data(eSet.brca.100)
#' @format expression set object
NULL


#' Example dataset
#' gene set projection input data example
#' 
#' @docType data
#' @name gspData
#' @usage data(gspData)
#' @format expression set object
#' @return \itemize{
#'  \item{"gsp.eSet"}{an eSet object}
#'  \item{"gsp.GeneSet"}{an S4 geneset object defined by \code{GeneSet-class}}
#' }
NULL

#' Example dataset
#' hyperenrichment input data example
#' 
#' @docType data
#' @name hyper
#' @usage data(hyper)
#' @format expression set object
#' @return \itemize{
#'  \item{"hyperGsets"}{an S4 geneset object defined by \code{GeneSet-class}}
#'  \item{"hyperSig"}{a list of character vectors (geneSets)}
#' }
NULL