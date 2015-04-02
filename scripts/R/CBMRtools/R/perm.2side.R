monotonize <- function( x, rnk=NULL )
{
  # ensure that the x values (sorted according to rnk) are
  # non-decreasing
  #
  x.srt <- if (is.null(rnk)) x else x[order(rnk)]
  for ( i in 2:length(x) ) x.srt[i] <- max(x.srt[i-1],x.srt[i])
  if (is.null(rnk)) x.srt else x.srt[rank(rnk)]
}
perm.2side.names <- c("p1","p2", "GC.p", "maxT", "fpr","fwer","fdr")

perm.2side.online <- function( obs, perm, dir=1, do.abs=F )
{
  # given an observed score vector and a permuted score vector,
  # calculate summary statistics for different types of p-values
  #
  # INPUT:
  #   obs - vector of observed scores (sorted from smallest to largest)
  #  perm - vector of permuted scores (matching the observed ones)
  #   dir - alternative being tested (1: "greater", -1: "less" )
  #
  i.p <- obs>=0
  i.n <- !i.p
  abs.prm <- abs(perm)
  abs.obs <- abs(obs)
    
  smry <- matrix( NA, length(obs), 7 )
  smry[,1]    <- as.numeric( dir*obs<=dir*perm )            # 1-sided p-value

  if ( do.abs ) {
    smry[,2] <- as.numeric( abs.obs<=abs.prm )              # 2-sided p-value (based on abs value)  
  }
  else {
    smry[i.p,2] <- as.numeric( obs[i.p]<=perm[i.p] )        # 2-sided p-value (based on min procedure)
    smry[i.n,2] <- as.numeric( obs[i.n]>=perm[i.n] )        # ..
  }
  perm.srt <- sort(perm)
  smry[i.p,3] <- as.numeric( obs[i.p]<=perm.srt[i.p] )      # 'rank-based' p-value
  smry[i.n,3] <- as.numeric( obs[i.n]>=perm.srt[i.n] )      # ..

  prm.srt <- monotonize( abs.prm, rnk=abs.obs )             # maxT p-value
  smry[,4] <- as.numeric(abs.obs<=prm.srt)                  #  ..
  
  smry[,5]    <- cumineq( abs.prm, obs=abs.obs, dir=2 )     # FPR
  smry[,6]    <- as.numeric( smry[,5]>0 )                   # FWER
  smry[,7]    <- smry[,5]                                   # FDR
  smry
}
perm.2side <- function( x, y, nperm=100, score, ngenes=NULL, seed=NULL,
                        alternative=c("two.sided","greater","less"), control=NULL,
                        balanced=FALSE, exhaustive=FALSE, equalized=FALSE, match.score=FALSE, do.abs=do.abs,
                        rnd=NULL, smoother=0, verbose=FALSE, debug=FALSE, ... )
{
  # Rank genes according to given score, and carry out permutation
  # test. 
  #
  # INPUT:
  #  - x           (m x n) matrix for m genes and n experiments
  #  - y           n-tuple of 'output' (can be class labels, or Surv object, etc)
  #  - nperm       number of permutation iterations
  #  - score       the score function to use to rank genes
  #  - ngenes      no. of genes to report
  #  - seed        random number generator seed
  #  - alternative hypothesis in the 1-sided p-value (greater: up in 0; less: up in 1)
  #  - balanced    computation based on balanced permutations (i.e., the
  #                permuted classes contain an equal proportion of the
  #                true classes)
  #  - equalized   set to true if want to carry out balanced permutation
  #                but the classes have unequal-size (not tested)
  #  - match.score if T, when balanced and unequal-sized classes, compute
  #                a matching observed score for each permuted score
  #  - ...         score-specific arguments
  #
  # OUPTUT
  #    (m x 9) matrix with columns
  #
  #      score (the observed score)
  #      p1 (1-sided p-value)
  #      etc.
  #
  alternative <- match.arg( alternative )
  score <- match.fun( score )
  #surv <- is.Surv(y)
  surv <- !is.null(dim(y))
  lev <- if ( surv ) c(0,1) else my.levels(y)
  ngenes <- min(ngenes,nrow(x))
  paired.score <- ("paired" %in% names(list(...)) && list(...)$paired )
  VERBOSE( verbose && paired.score, "  Paired analysis\n" )

  if ( ngenes<2 )
    stop( "ngenes must be >=2" )
  if ( surv  && nrow(y)!=ncol(x) )
    stop( "y (surv) must be same length as ncol(x)" )
  if ( !surv && length(y)!=ncol(x) )
    stop( "y must be same length as ncol(x)" )

  n <- if ( is.null(dim(y)) ) length(y) else nrow(y)

  VERBOSE( verbose, "  Computing observed scores .. " )
  x.obs <- score( x, y, ... )
  x.idx <- c(1:round(ngenes/2), (nrow(x)-(ngenes-round(ngenes/2))+1):nrow(x))
  VERBOSE( verbose, "sorting them .. " )
  x.idx <- order(x.obs)[x.idx]
  x.obs <- x.obs[x.idx]
  VERBOSE( verbose, "done.\n" )

  if ( nperm<=0 ) {
    return( x.obs )
  }  
  x.prm <- NULL
  x.prm <- matrix( 0, length(x.idx), length(perm.2side.names) )
  colnames(x.prm) <- perm.2side.names
  rownames(x.prm) <- names(x.obs)

  VERBOSE( verbose, "  Permutation test (", nperm, " iterations", sep="" )
  VERBOSE( balanced && verbose,  ", balanced" )
  VERBOSE( equalized && verbose, ", equalized" )
  VERBOSE( verbose, ") .. " )

  if ( !is.null(seed) )
    set.seed(seed)  
  cls.prm <- if ( surv )
    permute.binarray( y[,2], nperm=nperm, balanced=balanced, equalized=equalized,
                      exhaustive=exhaustive, control=control, verbose=verbose )
  else if ( !paired.score )
    permute.binarray( y, nperm=nperm, balanced=balanced, equalized=equalized,
                      exhaustive=exhaustive, control=control, verbose=verbose )
  else
    permute.paired(y, nperm=nperm, exhaustive=exhaustive, balanced=balanced, control=control, verbose=verbose)

  if ( nperm!=nrow(cls.prm) )
  {
    VERBOSE(verbose, paste("\n\tupdated number of permutations: ",
                           nrow(cls.prm)," (was ",nperm,")\n",sep=""))
    nperm <- nrow(cls.prm)
  }
  OBS <- x.obs
  percent <- pctstep <- max( .1, round(1/nperm,2) )

  dir <- if (alternative!="less") 1 else -1
  for ( i in 1:nperm )
  {
    y.prm   <- cls.prm[i,]
    idx.prm <- !is.na(y.prm)  # needed when balanced==TRUE and equalized==TRUE
    y.prm   <- y.prm[idx.prm] # ..
    
    if ( surv ) {             # managing a survival-type labeling
      Y.prm <- y[idx.prm]
      Y.prm[y.prm==1,2] <- 1; Y.prm[y.prm==1,1] <- y[idx.prm,1][y.prm==1]
      Y.prm[y.prm==0,2] <- 0; Y.prm[y.prm==0,1] <- y[idx.prm,1][y.prm==0]
      y.prm <- Y.prm
    }
    PRM <- if (paired.score)
      score(x[,y.prm], y[y.prm], ...)[x.idx]
    else
      score(x[,idx.prm], y.prm, ...)[x.idx]

    tmp <- perm.2side.online(OBS,PRM,dir=dir)
    if ( !all(dim(tmp)==dim(x.prm)) ) stop("incompatible dimensions for summary")
    x.prm <- x.prm + tmp

    if ( verbose & i>=nperm*percent ) {
      VERBOSE( verbose, percent*100,"% ", sep="" )
      percent <- round(percent + pctstep,1)
    }
  }
  VERBOSE( verbose, "\n" )

  ## complete computation of some p-values (see perm.2side.summary 
  ## ..for the "offline" computation of these statistics
  list( pval=p2ss.add( x.prm, x.obs, nperm, smoother=smoother ), nperm=nperm )
}
p2ss.add <- function( x.prm, x.obs, nperm, do.abs=FALSE, smoother=0 )
{
  # complete computation of some p-values (see perm.2side.summary for
  # ..the "offline" computation of these statistics
  #
  if (sum(p1.idx <- colnames(x.prm)=="p1")==0)   stop("missing 'p1' column'")
  if (sum(p2.idx <- colnames(x.prm)=="p2")==0)   stop("missing 'p2' column'")
  if (sum(gp.idx <- colnames(x.prm)=="GC.p")==0) stop("missing 'GC.p' column'")
  if (sum(mt.idx <- colnames(x.prm)=="maxT")==0) stop("missing 'maxT' column'")
  if (sum(fp.idx <- colnames(x.prm)=="fpr")==0)  stop("missing 'fpr' column'")
  if (sum(fd.idx <- colnames(x.prm)=="fdr")==0)  stop("missing 'fdr' column'")
  p.idx <- p1.idx | p2.idx
  
  x.prm[,!p.idx] <- x.prm[,!p.idx]/nperm
  x.prm[,p1.idx] <- x.prm[,p1.idx] + smoother
  x.prm[,p2.idx] <- if (do.abs)
    (x.prm[,p2.idx] + smoother)
  else
    2*(apply(cbind(x.prm[,p2.idx],nperm-x.prm[,p2.idx]),1,min)+smoother)
  x.prm[,p.idx]  <- x.prm[,p.idx]/(smoother+nperm)
  
  x.prm[,mt.idx] <- monotonize(x.prm[,mt.idx],rnk=-abs(x.obs))
  x.prm[,fp.idx] <- x.prm[,fp.idx]/nrow(x.prm)
  x.prm[,fd.idx] <- x.prm[,fd.idx]/cumineq(abs(x.obs), obs=abs(x.obs), dir=2)
  x.prm[,fd.idx] <- apply( cbind(x.prm[,fd.idx],rep(1,nrow(x.prm))),1,min)
  
  x.prm <- cbind(x.prm,fdr.test=x.prm[,fd.idx])
  x.prm[,fd.idx] <- pval2fdr(x.prm[,p2.idx])    
  x.prm <- cbind(score=x.obs,x.prm)
  
  x.prm <- x.prm[,match(c("score", "p1","p2","fdr","maxT"),
                        colnames(x.prm))]
  as.data.frame(x.prm,check.names=FALSE,stringsAsFactors=FALSE)
}
res.generate <- function(cls,seed=NULL,delta=10,sd=1,ngenes=500,do.plot=FALSE,pch=".")
{
  if (!is.null(seed)) set.seed(seed)
  if ( is.null(levels(cls))) levels(cls) <- sort(unique(cls))
  if ( length(levels(cls))!=2 ) stop( "cls must be binary" )
  n1 <- sum(cls==levels(cls)[1])
  n2 <- sum(cls==levels(cls)[2])

  deltas <- rev(seq(0,delta,delta/(ngenes-1)))
  ngenes <- length(deltas)
  x <- matrix(rnorm(ngenes*n1,mean=deltas,sd=sd),ngenes,n1)
  y <- matrix(rnorm(ngenes*n1,mean=deltas,sd=sd),ngenes,n1)
  out <- rbind(cbind(x, matrix(rnorm(ngenes*n1,sd=sd),ngenes,n1)),
               cbind(matrix(rnorm(ngenes*n1,sd=sd),ngenes,n1),y))
  colnames(out) <- c(paste("expt.1.",1:n1,sep=""),paste("expt.2.",1:n2,sep=""))
  rownames(out) <- c(paste("gene.1.",1:ngenes,sep=""),
                     paste("gene.2.",1:ngenes,sep=""))
  if (do.plot)
  {
    plot(apply(out[1:ngenes,1:n1],1,mean),
         apply(out[(ngenes+1):nrow(out),(n1+1):ncol(out)],1,mean), pch=pch, 
         xlab=paste("mean(class=",levels(cls)[1],")",sep=""),
         ylab=paste("mean(class=",levels(cls)[2],")",sep=""))
    points(apply(out[1:ngenes,(n1+1):ncol(out)],1,mean),
           apply(out[1:ngenes,1:n1],1,mean),
           pch=pch,col="red")
    points(apply(out[(ngenes+1):nrow(out),(n1+1):ncol(out)],1,mean),
           apply(out[(ngenes+1):nrow(out),1:n1],1,mean),
           pch=pch,col="green")
    abline(0,1)
  }
  out
}
# testing
#
if ( F )
{
  source( "~/dvlp/R/mysource.R")
  source( "~/dvlp/R/cumineq.R")
  source( "~/dvlp/R/gene.cluster.R")
  source( "~/dvlp/R/permute.array.R")
  source( "~/dvlp/R/perm.2side.R")
  
  x <- res.generate(nrow=50,ncol=20,delta=1, seed=111)
  heatmap(x$res,Rowv=NA,Colv=NA)
  
  y <- perm.2side( x$res, x$cls, nperm=1000, verbose=TRUE, score=snr, seed=111)
  source( "~/dvlp/R/perm.2side.old.R")
  y1 <- perm.2side( x$res, x$cls, nperm=1000, verbose=TRUE, score=snr, seed=111)
  all(y[,-6]==y1[,c(1,2,3,7,6,5,4,8)])
}
