monotonize <- function( x, rnk=NULL )
{
  # ensure that the x values (sorted according to rnk) are
  # non-decreasing
  #
  x.srt <- if (is.null(rnk)) x else x[order(rnk)]
  for ( i in 2:length(x) ) x.srt[i] <- max(x.srt[i-1],x.srt[i])
  if (is.null(rnk)) x.srt else x.srt[rank(rnk)]
}
perm.1side.names <- c("p", "GC.p", "maxT", "fpr", "fwer", "fdr" )

perm.1side.online <- function( obs, perm )
{
  ## given an observed score vector and a permuted score vector,
  ## calculate summary statistics for different types of p-values
  ##
  perm.srt <- sort(perm)
  abs.ord <- order(abs(obs))
  prm.monotonized <- -monotonize(-perm,length(perm):1)
    
  smry <- matrix( NA, length(obs), length(perm.1side.names) )
  smry[,1] <- as.numeric( obs>=perm )            # nominal p-value
  smry[,2] <- as.numeric( obs>=perm.srt )        # rank-based p-value
  smry[,3] <- as.numeric( obs>=prm.monotonized ) # maxT p-value
  smry[,4] <- cumineq( perm, obs=obs, dir=1 )    # FPR
  smry[,5] <- as.numeric( smry[,4]>0 )           # FWER
  smry[,6] <- smry[,4]                           # FDR
  smry
}
perm.1side <- function( x, y, score, nperm=100, ngenes=NULL, seed=NULL, smooth=0,
                        online=TRUE, control=NULL, rnd=NULL, debug=FALSE, verbose=FALSE, exhaustive=FALSE, ... )
{
  # rank genes according to given score, and carry out permutation test
  #
  # INPUT:
  #  - x       (m x n) matrix for m genes and n experiments
  #  - y       n-tuple of 'output' (can be class labels, or Surv object, etc)
  #  - nperm   number of permutation iterations
  #  - score   the score function to use to rank genes
  #  - ...     score-specific arguments
  #
  # OUPTUT
  #  (m x (nperm+2)) matrix with columns
  #
  #    perm_1 perm_2 .. perm_nperm order score
  #
  #  where
  #    perm_i is the permuted cox scores for the i-th iteration
  #    order  is the rank index of the genes according to their
  #           observed score
  #    score  is the sorted list of observed scores

  score <- match.fun( score )
  surv <- is.Surv(y)
  lev <- if ( surv ) c(0,1) else my.levels(y)
  ngenes <- min(ngenes,nrow(x))

  if ( ngenes<1 )
    stop( "ngenes must be >=1" )
  #if ( length(lev)>2 )
  #  stop( "y must be binary" )
  if ( surv  && nrow(y)!=ncol(x) )
    stop( "y must be same length as ncol(x)" )
  if ( !surv && length(y)!=ncol(x) )
    stop( "y must be same length as ncol(x)" )
  
  n <- if ( is.null(dim(y)) ) length(y) else nrow(y)
  
  VERBOSE( verbose, "Computing observed score .. " )
  x.obs <- score( x, y, verbose=verbose, ... )
  x.idx <- c(1:ngenes)
  x.idx <- order(x.obs)[x.idx]
  x.obs <- x.obs[x.idx]
  VERBOSE( verbose, "done.\n" )
  
  x.prm <- NULL
  x.prm <- matrix( 0, length(x.idx), length(perm.1side.names), dimnames=list(rownames(x)[x.idx],NULL) )
  colnames(x.prm) <- perm.1side.names

  #rownames(x.prm) <- names(x.obs)
  if ( !is.null(seed) ) set.seed(seed)
  VERBOSE( verbose, "Generating permuted class templates .. " )
  cls.prm <- if ( surv )
    permute.binarray( y[,2], nperm=nperm, balanced=F, equalized=F, control=control,
                      exhaustive=exhaustive, verbose=verbose )
  else
    permute.binarray( y, nperm=nperm, balanced=F, equalized=F, control=control,
                      exhaustive=exhaustive, verbose=verbose )

  if ( nperm!=nrow(cls.prm) )
  {
    VERBOSE(verbose, paste("updated number of permutations: ",
                           nrow(cls.prm)," (was ",nperm,")\n",sep=""))
    nperm <- nrow(cls.prm)
  }
  VERBOSE( verbose, "done.\n" )

  VERBOSE( verbose, "Permutation test (", nperm, " iterations) .. ", sep="" )
  percent <- pctstep <- max( .1, round(1/nperm,2) )
  for ( i in 1:nperm )
  {
    y.prm   <- cls.prm[i,]
    idx.prm <- !is.na(y.prm)  # needed when balanced==T and equalized==T
    y.prm   <- y.prm[idx.prm] # ..
    
    if ( surv ) {             # managing a survival-type labeling
      Y.prm <- y[idx.prm]
      Y.prm[y.prm==1,2] <- 1; Y.prm[y.prm==1,1] <- y[idx.prm,1][y.prm==1]
      Y.prm[y.prm==0,2] <- 0; Y.prm[y.prm==0,1] <- y[idx.prm,1][y.prm==0]
      y.prm <- Y.prm
    }
    PRM <- score(x[,idx.prm], y.prm, ...)[x.idx]

    tmp <- perm.1side.online(x.obs,PRM)
    if ( !all(dim(tmp)==dim(x.prm)) ) stop("incompatible dimensions for summary")
    x.prm <- x.prm + tmp
    
    if ( verbose & i>=nperm*percent ) {
      VERBOSE( verbose, percent*100,"% ", sep="" )
      percent <- round(percent + pctstep,1)
    }
  }
  VERBOSE( verbose, "done.\n" )

  list( pval=p1ss.add( x.prm, x.obs, nperm, smoother=smooth), nperm=nperm )
}
p1ss.add <- function( x.prm, x.obs, nperm, smoother=0 )
{
  p.idx  <- match.index('p',colnames(x.prm))
  fp.idx <- match.index('fpr',colnames(x.prm))
  fd.idx <- match.index('fdr',colnames(x.prm))
  mt.idx <- match.index('maxT',colnames(x.prm))
  
  x.prm[,-p.idx] <- x.prm[,-p.idx]/nperm
  x.prm[,p.idx]  <- x.prm[,p.idx] + smoother
  x.prm[,p.idx]  <- x.prm[,p.idx]/(smoother+nperm)
  x.prm[,mt.idx] <- monotonize(x.prm[,mt.idx])
  
  x.prm[,fd.idx] <- x.prm[,fp.idx]/cumineq(abs(x.obs), obs=abs(x.obs), dir=1)
  x.prm[,fd.idx] <- apply( cbind(x.prm[,fd.idx],rep(1,nrow(x.prm))),1,min)
  x.prm[,fp.idx] <- x.prm[,fp.idx]/nrow(x.prm)
  
  x.prm <- cbind(x.prm,fdr.test=x.prm[,fd.idx])
  x.prm[,fd.idx] <- pval2fdr(x.prm[,p.idx])
  rnames <- rownames(x.prm)
  x.prm <- cbind(score=x.obs,x.prm)
  rownames(x.prm) <- rnames

  x.prm <- x.prm[,match(c("score","p","fdr","maxT","fwer","fpr"),colnames(x.prm))]
  as.data.frame(x.prm,check.names=FALSE,stringsAsFactors=FALSE)
}
# source("~/dvlp/R/perm.1side.new.R")
# testing
#
if ( F )
{
  source( "~/dvlp/R/mysource.R")
  source( "~/dvlp/R/cumineq.R")
  source( "~/dvlp/R/gene.cluster.R")
  source( "~/dvlp/R/permute.array.R")
  source( "~/dvlp/R/perm.1side.R")
  
  x <- res.generate(nrow=50,ncol=20,delta=1, seed=111)
  heatmap(x$res,Rowv=NA,Colv=NA)
  
  y <- perm.1side( x$res, x$cls, nperm=1000, verbose=T, score=snr, online=T, seed=111)
  source( "~/dvlp/R/perm.1side.old.R")
  y1 <- perm.1side( x$res, x$cls, nperm=1000, verbose=T, score=snr, online=T, seed=111)
  all(y[,-4]==y1[,c(1,2,6,5,4,7,3)])
}
