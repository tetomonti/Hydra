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
##    name1 [1], name2 [1], name3 [2], ...
##  
##  [1] Bioinformatics Program, Boston University
##  [2] Section of Computational Biomedicine, Boston University  
#####################################################################################

#####################################################################################
## BEGIN documentation support (what follows are keyworded entries from which documentation pages will be extracted automatically)

#' replicateFiltering
#' 
#' \code{replicateFiltering} ... description ...
#'
#' @param <PARAM_1> parameter description
#' @param <PARAM_2> parameter description
#' ...
#'
#' @examples
#'
#' # comment
#' <CODE SNIPPET>
#' ...
#' @export

## END documentation support
#####################################################################################

## REPLICATE FILTERING
##
replicateFiltering <- function(
    dat,
    cls,
    nperm=0,
    control=NULL,
    seed=NULL,
    robust=(nperm==0),
    exhaustive=F,
    smooth=0,
    rnd=NULL,
    verbose=FALSE,
    skip.singleton=TRUE,
    min.sd=0)
{
  if ( !is.factor(cls) ) cls <- as.factor(cls)

  # making sure that labels start at 0
  #
  if ( any((match(cls,my.levels(cls))-1)!=cls) )
    stop("cls must be contiguous labels starting at 0")
  if ( is.null(levels(cls)) )
    stop( "levels(cls) cannot be NULL" )

  x <- dat@signal
  
  # eliminating singletons (i.e., samples w/o replicates), if any
  #
  if ( min(my.tabulate(cls))<2 )
  {
    if ( skip.singleton )
    {
      idx <- mmatch( which(my.tabulate(cls)<2)-1, cls)
      x <- x[,-idx]
      cls <- cls[-idx]
      cls <- match( cls, unique(cls) ) - 1
      levels(cls) <- unique(cls)

      if (!is.null(control)) {
        control <- control[-idx]
        control <- match( control, unique(control) ) - 1
        levels(control) <- unique(control)
      }
      VERBOSE( verbose, "skipped singleton samples {", paste( idx,collapse=", "), "}\n" )
    }
    else {
      stop("need at least two replicates per class")
    }
  }
  if ( nperm==0 )
  {
    VERBOSE( verbose, "Computing observed score " )
    VERBOSE( verbose && robust, "(median-based, slower) " )
    VERBOSE( verbose, ".. " )
    f.obs <- rgf.fstat( x, cls, robust=robust, min.sd=min.sd )
    VERBOSE( verbose, "done.\n" )
    if ( !is.null(rnd) )
      f.obs <- round(f.obs,rnd)

    df1 <- nlevels(cls)-1         # between-group degrees of freedom
    df2 <- sum(tabulate(cls)-1)   # within-group degrees of freedom
    pvl <- sapply(f.obs[,"Fstar"],pf,df1=df1,df2=df2,lower.tail=F)
    fdr <- pval2fdr(pvl,na.rm=T)
    f.obs <- cbind(score=f.obs[,"Fstar",drop=F],asymp.p=signif(pvl,2),asymp.fdr=signif(fdr,2))
    colnames(f.obs)[1] <- "score"

    dsc.idx <- match(rownames(f.obs),genenames(dat))
    if ( any(is.na(dsc.idx)) ) {
      warning("some probe names do not have matching descriptions")
    }
    else {
      f.obs <- data.frame(f.obs,description=dat@description[dsc.idx],stringsAsFactors=F)
    }    
    VERBOSE( verbose, "done.\n" )
    RGF <- list(nperm=0,pval=f.obs)
    return (RGF)
  }
  # ELSE ...
  #
  score <- function(x, cls, verbose=F){ -(rgf.fstat(x,cls,robust=robust,min.sd=min.sd)[,1]) }
  f.prm <- perm.1side( x, y=cls, score=score, nperm=nperm, online=online, seed=seed, smooth=smooth,
                       rnd=rnd, control=control, exhaustive=exhaustive, verbose=verbose )
  f.prm$pval <-
    f.prm$pval[,match(c("score","p","fdr","maxT","fwer"),colnames(f.prm$pval))]
  
  f.prm$pval[,1] <- -f.prm$pval[,1]
  
  if (!is.null(rnd) ) {
    f.prm$pval[,1] <- round(f.prm$pval[,1],rnd)
    f.prm$pval[,-1] <- signif(f.prm$pval[,-1],2)
  }
  dsc.idx <- match(rownames(f.prm$pval),genenames(dat))
  if ( any(is.na(dsc.idx)) ) {
    warning("some probe names are do not have matching descriptions")
  }
  else {
    f.prm$pval <- data.frame(f.prm$pval,description=dat@description[dsc.idx],stringsAsFactors=F)
  }
  f.prm    
}
## RGF FSTAT
##
rgfFstat <- function( X, cls, robust=F, min.sd=0 )
{
  # INPUT:
  #
  #       X - (m x n) matrix (m genes, n samples)
  #     cls - (n x 1) vector of class labels
  #  robust - use mean (default) or median (not implemented yet)
  
  #if ( any((match(cls,my.levels(cls))-1)!=cls) )
  #  stop("cls must be contiguous labels starting at 0")

  if ( !is.factor(cls) ) cls <- as.factor(cls)  

  levs <- my.levels(cls)
  nlev <- length(levs)
  tlev <- my.tabulate(cls)
  N <- length(cls)
  In <- rep(1,N)

  if ( robust )
  {
    MN <- apply( X, 1, median, na.rm=T )

    if ( max(tabulate(cls))>2 )
      mn <- sapply( levs, function(z) apply( X[,cls==z,drop=F], 1, median, na.rm=T ) )
    else {
      CLS <- sapply( levs, function(z){as.numeric(cls==z)} )
      mn <- t( t(X %*% CLS)/tlev )
    }
  }
  else
  {
    # overall mean
    #
    MN <- drop( (X %*% In) / N )

    # within-group means
    #
    CLS <- sapply( levs, function(z){as.numeric(cls==z)} )
    mn <- t( t(X %*% CLS)/tlev )
  }
  mnn <- matrix( NA, nrow(mn), N )
  for ( i in 1:nlev ) mnn[,cls==levs[i]] <- mn[,i]
  
  wt <- (X - mnn)^2%*%In / (N-nlev)      # MSE
  bt <- ((mn-MN)^2 %*% tlev) / (nlev-1)  # MSTR

  if (min.sd>0) {
    bt[bt<min.sd] <- min.sd
    wt[wt<min.sd] <- min.sd
  }
  Fstar <- bt/wt
  Fstar[bt==0 & wt==0] <- 0

  fstats <- data.frame(Fstar=Fstar[,1], B=bt[,1], W=wt[,1])
  return( fstats )
}
