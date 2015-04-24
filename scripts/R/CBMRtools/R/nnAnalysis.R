##  Copyright (c) 2013, 2014, Boston University. All rights reserved.
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
##    Stfeano Monti [1,2]
##
##  [1] Bioinformatics Program, Boston University
##  [2] Section of Computational Biomedicine, Boston University

## BEGIN documentation support (what follows are keyworded entries from which documentation pages will be extracted automatically)
##
#' nnAnalysis
#'
#' \code{nnAnalysis} rank rows (genes) by their 'similarity' to a
#' prototype profile and carries out permutation-based p-value
#' calculation
#'
#' @param dat NxM genes-by-sample matrix ('ExpressionSet' or 'resdata' or 'gctdata' object)
#' @param y M-tuple vector to use as the prototype profile
#' @param y.i alternatively, specify the row index in dat where to find the prototype profile
#' @param score similarity measure, must be one of 'pearson', 'spearman' or 'euclid'
#' @param nperm number of permutation iterations (default 0)
#' @param control constrained permutations (control for confounding factor)
#' @param seed random number generator seed (for reproducible results)
#' @param verbose verbose output
#' @param smooth smoother to prevent 0 p-values (default 1)
#'
#' @return list object with attributes 'pval' and 'nperm' corresponding to nn analysis permutation results
#' @examples
#'
#' require(Biobase)
#' data(eSet.brca.100)
#'
#'  ## passing the actual vector as argument
#'  ##
#'  NNout1 <- nnAnalysis(dat=eSet.brca.100,y=exprs(eSet.brca.100)[1,],score='pearson',
#'                       nperm=1000,seed=123,verbose=TRUE,smooth=1)
#'  head(NNout1$pval)
#'
#'  ## passing the index (within the matrix) as argument
#'  ##
#'  NNout2 <- nnAnalysis(dat=eSet.brca.100,y.i=1,score='pearson',nperm=1000,
#'                       seed=123,verbose=TRUE,smooth=1)
#'  head(NNout2$pval)
#'  all.equal(NNout1$pval,NNout2$pval) # TRUE
#' 
#'  ## controlling for a confounder (only discrete)
#'  ##
#'  CTL <- pData(eSet.brca.100)[,'ER_status']
#'  NNout3 <- nnAnalysis(dat=eSet.brca.100,y.i=1,score='pearson',nperm=1000,
#'                       seed=123,control=CTL,verbose=TRUE,smooth=1)
#'  head(NNout3$pval)
#'
#' @export

## END documentation support

########################################################################
#                              MAIN                                    #
########################################################################
##
## rank rows (genes) by their 'similarity' to a prototype profile
##
nnAnalysis <- function
(
 dat,          # NxM genes-by-sample matrix ('ExpressionSet' or 'resdata' or 'gctdata' object)
 y=NULL,       # M-tuple vector to use as the prototype profile
 y.i=NULL,     # alternatively, specify the row index in dat where to find the prototype profile
 score=c("pearson","spearman","euclid"),
               # 'similarity' measure
 nperm=0,      # number of permutation iterations
 control=NULL, # constrained permutations (control for confounding factor)
 seed=NULL,
 verbose=F,
 smooth=1      # smoother to prevent 0 p-values
 )
{
  ## handling of ExpressionSet data type
  ##
  if ( class(dat)=='ExpressionSet' )
  {
    dsc <- { # attempting to identify the description column
      if ( length(idx <- grep('symbol',colnames(fData(dat)),ignore.case=TRUE))==1 )
        sapply(fData(dat)[,idx],as.character)
      else if ( ncol(featureData(dat))>0 )
        sapply(fData(dat)[,1],as.character)
      else
        featureNames(dat)
    }      
    dat <- new('gctdata',signal=exprs(dat),description=dsc)
  }  
  # from here on, a gctdata/resdata object expected (see broad.file.formats.R)
  #
  score <- match.arg(score)
  SCORE <- switch(score,
                  euclid=function(dat,y) 
                  {
                    ## taking the negative, since it's a similarity measure
                    -sqrt(drop(t(t(dat)-y)^2 %*% rep(1,ncol(dat))))
                  },
                  spearman=function(dat,y)
                  {
                    ## map to ranks ..
                    dat <- t(apply(dat,1,rank))
                    y <- rank(y)

                    ## ..then pearson
                    dat <- ( dat-drop(fast.mean(dat)) )
                    y <- y - mean(y)
                    In <- rep(1,ncol(dat))
                    
                    drop(dat %*% y) / sqrt(drop(dat^2 %*% In) * (y^2 %*% In))
                  },
                  pearson=function(dat,y)
                  {
                    ## same as above (minus the rank mapping)
                    dat <- ( dat-drop(fast.mean(dat)) )
                    y <- y - mean(y)
                    In <- rep(1,ncol(dat))
                    
                    drop(dat %*% y) / sqrt(drop(dat^2 %*% In) * (y^2 %*% In))
                  },
                  stop("unrecognized score:",score))

  if ( !is.function(SCORE) ) stop( "undefined score" )
  if ( !xor(is.null(y),is.null(y.i)) ) stop( "either y or y.i must be specified" )
  if ( is.null(y) && (y.i>nrow(dat) || y.i<1) )  stop( "y.i out of range" )
  if ( is.null(y) )
  {
    y <- getSignal(dat)[y.i,]
    #dat <- dat[-y.i,]
  }
  if ( ncol(dat)!=length(y) ) stop( "*** ncol(dat)!=length(y) ***" )

  ## eliminating entries with no variation
  ##
  SD <- drop(fast.sd(getSignal(dat)))
  if ( any(SD<.00001) ) {
    SDrm <- SD<.00001
    dat <- subset.data(dat,genenames=genenames(dat)[!SDrm])
    warning('removed ', sum(SDrm), ' entries with no variation\n')
  }
  VERBOSE( verbose, "NN analysis based on ", score, "\n", sep="")
  
  nn.perm <- {
    if ( score=="euclid" )
      perm.1side(getSignal(dat), y, score=SCORE, nperm=nperm, seed=seed,
                 control=control, verbose=verbose, smoother=smooth )
    else
      perm.2side(getSignal(dat), y, score=SCORE, nperm=nperm, seed=seed,
                 control=control, verbose=verbose, smoother=smooth )
  }
  if ( nperm==0 )
    return(rev(nn.perm$pval))

  nn.perm$pval <- nn.perm$pval[rev(1:nrow(nn.perm$pval)),,drop=FALSE]
  
  list(pval=data.frame(nn.perm$pval,
           description=getDescription(dat)[match(rownames(nn.perm$pval),genenames(dat))],check.names=FALSE),
       nperm=nn.perm$nperm)
}
nnAnalysis.4gp <- function( gp.filename, gene.name, out.filename,
                             score=c("pearson","spearman","euclid"),
                             nperm=100, verbose=T )
{
  score <- match.arg(score)
  nperm <- as.numeric(nperm)
  prob <- as.numeric(prob)

  res <- read.res( gp.filename, verbose=verbose )
  x <- nnAnalysis( res, y.i=match(gene.name,rownames(res)),
                    score=score, nperm=nperm, verbose=verbose )

  my.write.matrix( x, justify="none", sep="\t", row.names=T, file=out.filename )
  VERBOSE( verbose, "output saved to file '", out.filename, "'.\n", sep="" )
  out.filename
}
## EXAMPLE OF USE (not assuming existence of CBMRtools package)
##
if ( FALSE )
{
  rm(list=ls())
  CBMGIT <- Sys.getenv('CBMGIT')
  if (CBMGIT=="") stop( "Use 'setenv CBMGIT ..' to set CBMgithub's base directory" )

  source(paste(CBMGIT,"scripts/R/CBMRtools/R/misc.R",sep="/"))
  source(paste(CBMGIT,"scripts/R/CBMRtools/R/misc.math.R",sep="/"))
  source(paste(CBMGIT,"scripts/R/CBMRtools/R/broad.file.formats.R",sep="/"))
  source(paste(CBMGIT,"scripts/R/CBMRtools/R/permute.array.R",sep="/"))
  source(paste(CBMGIT,"scripts/R/CBMRtools/R/perm.1side.R",sep="/"))
  source(paste(CBMGIT,"scripts/R/CBMRtools/R/perm.2side.R",sep="/"))
  source(paste(CBMGIT,"scripts/R/CBMRtools/R/diffanal.scores.R",sep="/"))
  source(paste(CBMGIT,"scripts/R/CBMRtools/R/nnAnalysis.R",sep="/"))
  
  require(Biobase)
  
  eSet.brca.100 <- load.var(paste(CBMGIT,"scripts/R/CBMRtools/data/eSet.brca.100.rda",sep='/'))

  ## passing the actual vector as argument
  ##
  NNout1 <- nnAnalysis(dat=eSet.brca.100,y=exprs(eSet.brca.100)[1,],score='pearson',nperm=1000,seed=123,verbose=TRUE,smooth=1)
  head(NNout1$pval)

  ## passing the index (within the matrix) as argument
  ##
  NNout2 <- nnAnalysis(dat=eSet.brca.100,y.i=1,score='pearson',nperm=1000,seed=123,verbose=TRUE,smooth=1)
  head(NNout2$pval)
  all.equal(NNout1$pval,NNout2$pval) # should be TRUE

  ## controlling for a confounder (only discrete)
  ##
  CTL <- pData(eSet.brca.100)[,'gender']
  NNout3 <- nnAnalysis(dat=eSet.brca.100,y.i=1,score='pearson',nperm=1000,seed=123,control=CTL,verbose=TRUE,smooth=1)
  head(NNout3$pval)
}
