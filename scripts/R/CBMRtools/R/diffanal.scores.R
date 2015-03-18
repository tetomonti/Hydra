#source( "~/dvlp/R/levels.R" )
#source( "~/dvlp/R/permute.array.R" )
#source( "~/dvlp/R/robust.sd.R" )
#source( "~/dvlp/R/plot.markers.R" )
#source( "~/dvlp/R/verbose.R" )

# function: T.SCORE
#
t.score <- function( x, cls=NULL, y=NULL, robust=FALSE, paired=FALSE, generalized=FALSE,
                     do.test=F, verbose=F, var.equal=F, min.sd=NULL, 
                     alternative=c("two.sided","greater","less") )
{
  # INPUT:
  #    x - m x n1 matrix (genes by experiments, condition 1)
  #    y - m x n2 matrix (genes by experiments, condition 2)
  #  OR
  #    x - m x n  matrix (genes by experiments, condition 1 & 2)
  #  cls - n vector of class labels
  #
  # OUTPUT:
  #  snr - m vector (positive are upregulated for x, or for
  #        lower label -- condition 1 -- when cls is specified)
  
  # some checks on the input
  #
  alternative <- match.arg(alternative)

  if ( is.null(y) & is.null(cls) )
    stop( "must specify either y or cls" )

  if ( is.null(y) )
  {
    lev <- sort(unique(cls))
    if ( ncol(x)!=length(cls) )
      stop( "ncol(x) must be same as length(cls)" )
    if ( length(lev)>2 )
      stop( "cls must be binary" )
    y <- x[,cls==lev[2],drop=F]
    x <- x[,cls==lev[1],drop=F]
  }  
  if ( nrow(x)!=nrow(y) ) stop( "x and y must be of same length\n" )
  if ( ncol(x)<4 ) warning( "x has less than 4 observations\n" )
  if ( ncol(y)<4 ) warning( "y has less than 4 observations\n" )

  score <- NULL
  
  # paired score
  #
  if ( paired )
  {
    if ( generalized )
    {
      if ( (ncol(x)>ncol(y) & ncol(x)%%ncol(y)) |
           (ncol(y)>ncol(x) & ncol(y)%%ncol(x)) )
        stop( "x and y must be have column numbers multiple of each other\n" )

      d <- NULL

      if ( ncol(x)>ncol(y) )
        for ( i in 0:(ncol(x)/ncol(y)-1) ) {
          idx <- (i*ncol(y)+1):((i+1)*ncol(y))
          VERBOSE(verbose, "comparing x[:",idx[1],":",idx[length(idx)],
                  "] to y[1:",ncol(y),"]","\n",sep="")
          d <- cbind( d, x[,idx,drop=F]-y )
        }
      else if ( ncol(y)>ncol(x) )
        for ( i in 0:(ncol(y)/ncol(x)-1) ) {
          idx <- (i*ncol(x)+1):((i+1)*ncol(x))
          VERBOSE(verbose, "comparing x[1:",ncol(x),"] to y[",
                  idx[1],":",idx[length(idx)],"]","\n",sep="")
          d <- cbind( d, x-y[,idx,drop=F] )
        }
      else
        d <- x-y
    }
    else
    {
      if ( ncol(x)!=ncol(y) )
        stop( "x and y must have same number of columns\n" )

      d <- x-y
    }
    VERBOSE( verbose, "computing paired t.score .." )
    if ( robust ) {
      stop( "robust paired t.score not implemented yet" )
    }
    else {
      if (do.test) {
        score <- cbind(score=drop((fast.mean(d)*sqrt(length(cls)))/fast.sd(d)),
                       p.value=drop(apply(d, 1, function(z) t.test(z)$p.value)))
        rownames(score) <- rownames(d)
      }
      else {
        score <- (drop(fast.mean(d))*sqrt(ncol(d)))/drop(fast.sd(d))
        names(score) <- rownames(d)
      }
    }
    VERBOSE( verbose, " done.\n" )
    return( score )
  }
  # ELSE not paired
  #
  x.idx <- 1:ncol(x)

  VERBOSE( verbose, ifelse( robust, "\tWilcoxon test .. ", "\tt test .. " ) )
  
  n1 <- (ncol(x)); if (n1<2) stop( "need at least 2 obs per class" )
  n2 <- (ncol(y)); if (n2<2) stop( "need at least 2 obs per class" )
  cls <- c( rep(1,n1), rep(0,n2) ); cls <- cbind( cls, 1-cls )
  x <- cbind(x,y)
  
  if ( robust )
  {
    rnk <- t(apply(-x,1,rank))
    score <-  (drop(rnk[,1:n1] %*% rep(1,n1)) - as.double(n1) * (as.double(n1)+1)/2)

    if (do.test) {
      warning( "wilcox.test p-value not implemented yet, ignoring" )
    }
  }
  ## FAST COMPUTATION BASED ON MATRIX MULTIPLICATION
  ##
  s  <- x %*% cls
  s2 <- x^2 %*% cls
  s2[,1] <- (s2[,1] - (s[,1]^2)/n1) / (n1-1) # variance in 1st class 
  s2[,2] <- (s2[,2] - (s[,2]^2)/n2) / (n2-1) # variance in 2nd class 
  s[,1] <- s[,1]/n1
  s[,2] <- s[,2]/n2
  
  if ( !is.null(min.sd) ) {
    min.var <- min.sd*min.sd
    s2[s2[,1]<min.var,1] <- min.var
    s2[s2[,2]<min.var,2] <- min.var
  }    
  stderr <- if (var.equal)
    sqrt( (((n1-1)*s2[,1] + (n2-1)*s2[,2])/(n1+n2-2)) * (1/n1+1/n2) )
  else
    sqrt( s2[,1]/n1 + s2[,2]/n2 )
  
  score <- (s[,1]-s[,2]) / stderr
  
  if ( do.test )
  {
    df <- if ( var.equal ) # degrees of freedom
      n1+n2-2
    else
      stderr^4 / ( (s2[,1]/n1)^2/(n1-1) + (s2[,2]/n2)^2/(n2-1)) # Welch approximation of df
    
    pval <- if (alternative == "less") {
      pt(score, df=df)
    }
    else if (alternative == "greater") {
      pt(score, df=df, lower.tail=F)
    }
    else {
      2 * pt(-abs(score), df=df)
    }
    score <- cbind( score=score, p.value=pval )
  }
  VERBOSE( verbose, "done.\n" )
  return( score )
}
# function: FIX SD
#
fix.sd <- function( s, m, s.percent=0.2 )
{
  ## function used in SNR computation, to threshold stdev as done in
  ## GeneCluster (for backward compatibility)
  ##
  if ( is.vector(s) || is.matrix(s) )
  {
    if ( ( is.vector(s) && length(s)!=length(m) ) ||
         ( is.matrix(s) && all(dim(s)!=dim(m)) ) )
      stop("s and m must be same length")

    abs.m <- abs( m )
    min.s <- abs.m * s.percent 
    min.s[min.s<s] <- s[min.s<s]
    min.s[min.s<=0] <- 0.1
    return( min.s )
  }
  else
  {
    abs.m <- abs( m )
    min.s <- s.percent * abs.m
    if ( min.s<s ) { min.s <- s }
    if ( min.s==0 ) { min.s <- 0.1 }
    return( min.s )
  }
}
## function: SNR
##
snr <- function( x, cls=NULL, y=NULL, robust=FALSE, fix=TRUE, gc=fix, min.sd=0.0,
                 paired=FALSE, generalized=FALSE, do.test=FALSE, verbose=FALSE )
{
  # INPUT:
  #    x - m x n1 matrix (genes by experiments, 1st condition)
  #    y - m x n2 matrix (genes by experiments, 2nd condition)
  #  OR
  #    x - m x n  matrix (genes by experiments, both conditions)
  #  cls - n vector of class labels
  #
  # OUTPUT:
  #  snr - m vector (positive are upregulated for x, or for
  #        lower label -- 1st condition -- when cls is specified)
  
  # some checks on the input
  #
  if ( is.null(y) && is.null(cls) )
    stop( "must specify either y or cls" )
  if ( fix && min.sd>0 )
    stop( "can't use both sd's fix and min.sd" )
  
  if ( is.null(y) )
  {
    lev <- sort(unique(cls))
    if ( ncol(x)!=length(cls) )
      stop( "ncol(x) must be same as length(cls)" )
    if ( length(lev)>2 )
      stop( "cls must be binary" )
    y <- x[,cls==lev[2],drop=FALSE]
    x <- x[,cls==lev[1],drop=FALSE]
  }  
  if ( nrow(x)!=nrow(y) ) stop( "x and y must be of same length\n" )
  if ( ncol(x)<4 ) warning( "x has less than 4 observations\n" )
  if ( ncol(y)<4 ) warning( "y has less than 4 observations\n" )

  # paired SNR
  #
  if ( paired )
  {
    if ( generalized )
    {
      if ( (ncol(x)>ncol(y) & ncol(x)%%ncol(y)) |
           (ncol(y)>ncol(x) & ncol(y)%%ncol(x)) )
        stop( "x and y must have column numbers multiple of each other\n" )

      d <- NULL

      if ( ncol(x)>ncol(y) )
        for ( i in 0:(ncol(x)/ncol(y)-1) ) {
          idx <- (i*ncol(y)+1):((i+1)*ncol(y))
          VERBOSE(verbose, "comparing x[:",idx[1],":",idx[length(idx)],
                  "] to y[1:",ncol(y),"]","\n",sep="")
          d <- cbind( d, x[,idx,drop=FALSE]-y )
        }
      else if ( ncol(y)>ncol(x) )
        for ( i in 0:(ncol(y)/ncol(x)-1) ) {
          idx <- (i*ncol(x)+1):((i+1)*ncol(x))
          VERBOSE(verbose, "comparing x[1:",ncol(x),"] to y[",
                  idx[1],":",idx[length(idx)],"]","\n",sep="")
          d <- cbind( d, x-y[,idx,drop=FALSE] )
        }
      else
        d <- x-y
    }
    else
    {
      if ( ncol(x)!=ncol(y) )
        stop( "x and y must have same number of columns\n" )

      d <- x-y
    }
    VERBOSE( verbose, "computing paired SNR .." )
    if ( robust ) {
      m <- apply( d, 1, median )
      s <- apply( d, 1, mad )
    }
    else {
      m <- apply( d, 1, mean )
      s <- apply( d, 1, sd )
    }
    if ( fix ) s <- fix.sd(s,m,s.percent=.05)
    snr <- m/s
    VERBOSE( verbose, " done.\n" )
    
    if ( do.test )
    {
      VERBOSE( verbose, "computing asymptotic p-values (n=", ncol(d), ") ..", sep="" )
      p.value <- apply( d,1,function(z){ t.test(z)$p.value } )
      return ( cbind(score=snr,p.value=p.value) )
      VERBOSE( verbose, " done.\n" )
    }
    return( snr )
  }
  # ELSE not paired
  #
  snr <- NULL
  if ( robust )
  {
    x.m <- apply( x, 1, median )
    x.s <- if ( gc ) apply( x, 1, sd) else apply( x, 1, mad )
    y.m <- apply( y, 1, median )
    y.s <- if ( gc ) apply( y, 1, sd ) else apply( y, 1, mad )

    # thresholding as implemented in GeneCluster (to handle zero variance)
    #
    if ( fix )
    {
      x.s <- fix.sd(x.s,x.m)
      y.s <- fix.sd(y.s,y.m)
    }
    if ( min.sd>0 )
    {
      x.s[x.s<min.sd] <- min.sd
      y.s[y.s<min.sd] <- min.sd
    }      
    snr <- (x.m-y.m)/(x.s+y.s)
  }
  else
  {
    n1 <- ncol(x)
    n2 <- ncol(y)
    cls <- c( rep(1,n1), rep(0,n2) ); cls <- cbind( cls, 1-cls )
    X <- cbind(x,y)
    
    s  <- X %*% cls
    s2 <- X^2 %*% cls
    s2[,1] <- (s2[,1] - (s[,1]^2)/n1) / (n1-1); s2[s2[,1]<0,1] <- 0
    s2[,2] <- (s2[,2] - (s[,2]^2)/n2) / (n2-1); s2[s2[,2]<0,2] <- 0
    s2[,1] <- sqrt( s2[,1] )
    s2[,2] <- sqrt( s2[,2] )
    s[,1] <- s[,1]/n1
    s[,2] <- s[,2]/n2
    
    # thresholding as implemented in GeneCluster (to handle zero variance)
    #
    if ( fix ) {
      s2[,1] <- fix.sd(s2[,1],s[,1])
      s2[,2] <- fix.sd(s2[,2],s[,2])
    }    
    if ( min.sd>0 )
    {
      s2[s2[,1]<min.sd,1] <- min.sd
      s2[s2[,2]<min.sd,2] <- min.sd
    }      
    snr <- (s[,1]-s[,2]) / (s2[,1]+s2[,2])
  }
  if ( do.test )
  {
    VERBOSE( verbose, "computing asymptotic p-values (n=[",
             paste(ncol(x),ncol(y),sep="x"), "]) ..", sep="" )
    p.value <- apply( cbind(x,y), 1,
                      function(z){t.test(z[1:ncol(x)],z[(ncol(x)+1):length(z)],var.equal=TRUE)$p.value} )
    return ( cbind(score=snr, p.value=p.value) )
    VERBOSE( verbose, " done.\n" )
  }
  else {
    return( snr )
  }
}
# function: CLS STATS
#
cls.stats <- function( x, cls=NULL, control=NULL, rnd=NULL, robust=FALSE, unlog=FALSE, gc=FALSE, paired=FALSE )
{
  if ( is.null(cls) )
    cls <- rep(0,ncol(x))
  lev <- sort(unique(cls))
  if (paired && length(lev)!=2) stop("binary cls expected w/ paired data")
  if ( is.null(levels(cls)) ) levels(cls) <- lev
  CLS <- sapply( 1:length(lev), function (z) as.numeric(cls==lev[z]) )
  nc  <- apply(CLS,2,sum)
  if (unlog) x <- 2^x

  if ( robust )
  {
    SCL <- if (gc) sd else mad
    s <- sapply( lev, function(z) apply(x[,cls==z,drop=FALSE],1,median) )
    s2 <- sapply( lev, function(z) apply(x[,cls==z,drop=FALSE],1,SCL) )
  }
  else
  {
    s  <- x %*% CLS
    s2 <- x^2 %*% CLS
    s2 <- sqrt( t((t(s2) - t(s^2)/nc) / (nc-1)) )
    s  <- t(t(s)/nc)
  }
  fold.chg <- {
    if ( paired ) {
      if (robust) 
        2^apply(log2(x[,cls==lev[1]])-log2(x[,cls==lev[2]]),1,median)
      else
        2^(drop(fast.mean(log2(x[,cls==lev[1]])-log2(x[,cls==lev[2]]))))
    }
    else {
      s[,1]/s[,2]
    }
  }
  if ( any(fold.chg<1) ) {
    fold.chg[fold.chg<1] <- 1/fold.chg[fold.chg<1]
  }
  ctr <- if (robust) "median" else "mean"
  scl <- if (robust) "mad" else "stdv"
  colnames(s) <- paste(ctr,levels(cls),sep=".")
  colnames(s2) <- paste(scl,levels(cls),sep=".")

  stats <- cbind(fold.chg=fold.chg,s,s2)
  if ( !is.null(rnd) ) {
    stats <- round(stats,rnd)
  }
  stats
}
