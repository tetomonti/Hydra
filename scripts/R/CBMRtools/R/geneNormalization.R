#####################################################################################
##  Copyright (c) 2013-2015, Boston University. All rights reserved.
##  
##  Redistribution and use in source and binary forms, with or without
##  modification, are permitted provided that the following conditions are met: 
##  
##  1. Redistributions of source code must retain the above copyright notice, this
##     list of conditions and the following disclaimer. 
##  2. Redistributions in binary form must reproduce the above copyright notice,
##     this list of conditions and the following disclaimer in the documentation
##     and/or other materials provided with the distribution. 
##  
##  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
##  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
##  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
##  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
##  ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
##  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
##  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
##  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
##  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
##  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##  
##  The views and conclusions contained in the software and documentation are those
##  of the authors and should not be interpreted as representing official policies, 
##  either expressed or implied, of Boston University.
##  
##  Authors:
##    Stefano Monti [1,2]
##  
##  [1] Bioinformatics Program, Boston University
##  [2] Section of Computational Biomedicine, Boston University  
#####################################################################################

#####################################################################################
## BEGIN documentation support

#' geneNormalization
#' 
#' \code{geneNormalization} gene-specific normalization, so that each
#' gene has the same mean/sd in two datasets
#'
#' @param norm dataset to normalize
#' @param base dataset with respect to which norm will be normalized
#' @param norm.sub subset of norm's samples to be used as reference
#' @param base.sub subset of base's samples to be used as reference 
#' @param robust use median/MAD in place of mean/sd
#'
#' @examples
#'
#' data(eSet.brca.100)
#' set.seed(123) # for reproducible results
#' dat.idx <- sort(sample(ncol(eSet.brca.100),ncol(eSet.brca.100)/2))
#' 
#' N <- eSet.brca.100[, dat.idx]
#' B <- eSet.brca.100[,-dat.idx]
#' 
#' set.seed(123)
#' N.sub <- sampleNames(N)[sort(sample(ncol(N),ncol(N)/2))]
#' B.sub <- sampleNames(B)[sort(sample(ncol(B),ncol(B)/2))]
#'
#' # normalize using all samples
#' NORMall <- geneNormalization(exprs(N),exprs(B))
#' plot(apply(exprs(B),1,mean),apply(NORMall,1,mean),xlab="base",ylab="norm",main="all")
#' all.equal(apply(exprs(B),1,mean),apply(NORMall,1,mean)) # TRUE
#'  
#' # normalize using a subset of samples in both datasets
#' NORMsub <- geneNormalization(exprs(N),exprs(B),norm.sub=N.sub,base.sub=B.sub)
#'
#' # the overall gene means are different
#' plot(apply(exprs(B),1,mean),apply(NORMsub,1,mean),
#'      xlab="base",ylab="norm",main="sub")
#' all.equal(apply(exprs(B),1,mean),apply(NORMsub,1,mean)) # FALSE
#'
#' # but they are the same within the subsets
#' plot(apply(exprs(B[,B.sub]),1,mean),apply(NORMsub[,N.sub],1,mean),
#'      xlab="base",ylab="norm",main="within sub")
#' all.equal(apply(exprs(B[,B.sub]),1,mean),apply(NORMsub[,N.sub],1,mean)) # TRUE
#'  
#' @export

## END documentation support
#####################################################################################
geneNormalization <- function( norm, base, norm.sub=NULL, base.sub=NULL, robust=FALSE )
{
  # normalize norm so that each gene has same "mean" and "stdev" as
  # the corresponding gene in base
  #
  # INPUT:
  #       norm - m x n1 matrix (m "genes", n1 "experiments")
  #       base - m x n2 matrix (m "genes", n2 "experiments")
  #   norm.sub - vector of sample names in norm to use for normalization
  #   base.sub - vector of sample names in base to use for normalization
  #
  # OUTPUT
  #  norm.norm - m x n1 matrix (m "genes", n1 "experiments")
  #
  if ( nrow(norm)!=nrow(base) ) {
    stop( "norm and base must have same number of rows" )
  }
  if ( !is.null(norm.sub) ) {
    if ( length(norm.sub)<2 )
      stop( "using less than 2 samples in norm subsetting" )
    if ( length(norm.sub)<3 )
      warning( "using less than 3 samples in norm subsetting" )
  }
  if ( !is.null(base.sub) ) {
    if ( length(base.sub)<2 )
      stop( "using less than 2 samples in base subsetting" )
    if ( length(base.sub)<3 )
      warning( "using less than 3 samples in base subsetting" )
  }
  if ( any(rownames(norm)!=rownames(base)) )
    stop( "norm and base must contain same genes in same order" )

  location <- if(robust) median else mean
  scale <- if(robust) mad else sd

  if ( is.null(norm.sub) ) {
    norm.sub <- colnames(norm)
  }
  if ( is.null(base.sub) ) {
    base.sub <- colnames(base)
  }
  if ( any(is.na(norm.idx <- match(norm.sub,colnames(norm)))) )
    stop( "samples in norm.sub missing from norm" )
  if ( any(is.na(base.idx <- match(base.sub,colnames(base)))) )
    stop( "samples in base.sub missing from base" )
  
  norm.mn <- apply( norm[,norm.idx], 1, location )
  norm.sd <- apply( norm[,norm.idx], 1, scale )
  base.mn <- apply( base[,base.idx], 1, location )
  base.sd <- apply( base[,base.idx], 1, scale )

  norm <- ( (norm - norm.mn) / norm.sd ) * base.sd + base.mn
  return( norm )
}
#######################################################################
# example of use (w/o pacakge installation)
#######################################################################
if ( F )
{
  CBMGIT <- Sys.getenv('CBMGIT')
  if (CBMGIT=="") stop( "Use 'setenv CBMGIT ..' to set CBMgithub's base directory" )
  source( paste(CBMGIT, "scripts/R/CBMRtools/R/geneNormalization.R", sep="/") )
  require(Biobase)
    
  load(paste(CBMGIT,'/scripts/R/CBMRtools/data/eSet.brca.100.rda',sep=''))
  set.seed(123) # for reproducible results
  dat.idx <- sort(sample(ncol(eSet.brca.100),ncol(eSet.brca.100)/2))
  
  N <- eSet.brca.100[, dat.idx]
  B <- eSet.brca.100[,-dat.idx]

  set.seed(123)
  N.sub <- sampleNames(N)[sort(sample(ncol(N),ncol(N)/2))]
  B.sub <- sampleNames(B)[sort(sample(ncol(B),ncol(B)/2))]

  ## normalize using all samples
  NORMall <- geneNormalization(exprs(N),exprs(B))
  plot(apply(exprs(B),1,mean),apply(NORMall,1,mean),xlab="base",ylab="norm",main="sub")
  all.equal(apply(exprs(B),1,mean),apply(NORMall,1,mean)) # TRUE
  
  ## normalize using a subset of samples in both datasets
  NORMsub <- geneNormalization(exprs(N),exprs(B),norm.sub=N.sub,base.sub=B.sub)
  plot(apply(exprs(B),1,mean),apply(NORMsub,1,mean),xlab="base",ylab="norm",main="all")
  all.equal(apply(exprs(B),1,mean),apply(NORMsub,1,mean)) # FALSE

  plot(apply(exprs(B[,B.sub]),1,mean),apply(NORMsub[,N.sub],1,mean),xlab="base",ylab="norm",main="all")
  all.equal(apply(exprs(B[,B.sub]),1,mean),apply(NORMsub[,N.sub],1,mean)) # TRUE
}
