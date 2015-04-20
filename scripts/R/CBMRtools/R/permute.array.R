# function JOIN CROSSPROD
#
join.crossprod <- function( mx1, mx2 )
{
  if ( !is.matrix(mx1) )
    mx1 <- t(as.matrix(mx1))
  if ( !is.matrix(mx2) )
    mx2 <- t(as.matrix(mx2))

  mx.out <- matrix( NA, nrow(mx1)*nrow(mx2), ncol(mx1)+ncol(mx2) )
  base <- 0
  for ( i in 1:nrow(mx1) ) {
    mx.out[(base+1):(base+nrow(mx2)),] <- cbind(matrix(mx1[i,],nrow(mx2),ncol(mx1),byrow=T),mx2)
    base <- base+nrow(mx2)
  }
  mx.out
}
# function PERMUTE PAIRED
#
permute.paired <- function( cls, nperm=1, balanced=F, exhaustive=F, control=NULL, verbose=F, seed=NULL )
{
  # return matrix whereby each row represents the indeces for a paired
  # rearrangement of columns. It assumes cls to have binary and
  # contiguous labels.  Example: assume cls=c(0,0,0,1,1,1). Then the
  # non-permuted index vector would be c(1,2,3,4,5,6), which says that
  # 1 is paired with 4, 2, with 5, etc. A plausible (4x6) matrix of
  # permuted indeces could be:
  #   1,5,3,4,2,6
  #   4,5,3,1,2,6
  #   4,2,6,1,5,3
  #   4,2,3,1,5,6
  #
  lev  <- my.levels(cls)
  m    <- length(cls)
  idx1 <- (1:m)[cls==lev[1]]
  idx2 <- (1:m)[cls==lev[2]]
  nc   <- length(idx1)

  # some checks on input
  #
  if ( length(lev)>2 )
    stop( "paired permutation requires binary class" )
  if ( length(idx1)!=length(idx2) )
    stop( "paired permutation requires balanced classes" )
  if ( any(idx1!=(1:nc)) & any(idx2!=(1:nc)) )
    stop( "paired samples must be contiguous" )
  if ( nc<2 )
    stop( "need at least two elements per class" )
  if ( nc==2 & nperm>2 )
    stop( "with two elements per class, only two permuations available" )


  if ( !is.null(seed) )
    set.seed(seed)

  # determining the exhaustive number of possible permutations
  #
  ntot <- 0
  k1 <- if (balanced) nc %/% 2 else 1
  k2 <- if (balanced) k1 + (nc %% 2) else nc
  for ( k in k1:k2 ) ntot <- ntot+choose(nc,k)
  if (!balanced) ntot <- ntot+1
  
  if ( exhaustive && ntot>nperm) {
    VERBOSE(T, "Number of exhaustive permutations is too large:", ntot, "(ignored)\n")
    exhaustive <- F
  }
  else if ( !exhaustive && nperm>ntot ) {
    VERBOSE(T, "Number of possible permutations less than needed:", ntot, "(reduced)\n")
    nperm <- ntot
    exhaustive <- T
  }
  # EXHAUSTIVE set of permutations
  #
  ord <- 1:m
  if ( exhaustive )
  {
    VERBOSE( verbose, "Number of permutations:", ntot, "\n" )
    
    lbl <- matrix( rep(ord,ntot), ntot, m, byrow=T )
    offset <- 0
    for ( k in k1:k2 )
    {
      idx1.perm <- t( combn(nc,k) )
      for ( i in (1:nrow(idx1.perm)) )
      {
        idx <- idx1.perm[i,]
        lbl[i+offset,idx] <- ord[idx+nc] # switch k indeces btw class1 ..
        lbl[i+offset,idx+nc] <- ord[idx] # ..and class2
      }
      offset <- offset + i
    }
    return(lbl)
  }
  # ELSE (non-exhaustive set of permuations) ..
  #
  lbl <- matrix( rep(1:m,nperm), nperm, m, byrow=T )
  
  if ( balanced )
  {
    k1 <- nc %/% 2
    k2 <- nc %% 2
    for ( i in (1:nperm) )
    {
      k <- k1 + sample(c(0,1),1)*k2
      idx1.perm <- sort(sample(idx1, k))
      idx2.perm <- idx1.perm + nc
      lbl[i,idx1.perm] <- ord[idx2.perm]
      lbl[i,idx2.perm] <- ord[idx1.perm]
    }
  }
  else # NOT balanced
  {
    for ( i in (1:nperm) )
    {
      k <- sample( 1:(nc-1), 1 ) 
      idx1.perm <- sort(sample( idx1, k ))
      idx2.perm <- idx1.perm + nc
      lbl[i,idx1.perm] <- ord[idx2.perm]
      lbl[i,idx2.perm] <- ord[idx1.perm]
    }
  }
  return( lbl )
}
permute.binarray <- function(cls, nperm=1, balanced=F, equalized=F, seed=NULL,
                             exhaustive=F, control=NULL, verbose=F )
{
  if ( !is.null(dim(cls)) )
    if (ncol(cls)==2)
      cls <- cls[,2]
    else
      stop( "multi-dimensional class label must have 2 columns" )
  lev <- my.levels(cls)
  m  <- length(cls)
  perm <- matrix( NA, nrow=nperm, ncol=length(cls) )
  if ( !is.null(seed) ) set.seed(seed)

  # non-binary response variable
  #
  if ( length(lev)>2 )
  {
    if ( is.null(control) ) {
      for ( i in (1:nperm) ) {
        perm[i,] <- sample(cls,m)
      }
      return(perm)
    }
    # control for confounder
    #
    else
    {
      control <- many2one(control)
      ctl.levs <- sort(unique(control))

      VERBOSE( verbose, "\tcontrolling for ", length(ctl.levs), "-level confounder\n", sep="" )

      idx <- lapply( ctl.levs, function(z) which(control==z) )
      jar <- lapply( idx, function(z) cls[z] )

      nperm.tot <- 0
      for ( i in 1:length(jar) ) {
        nperm.tot <- nperm.tot + (tmp <- lmchoose(my.tabulate(jar[[i]])))
        VERBOSE( verbose, "  nchoose[",ctl.levs[i],"]:   ", exp(tmp), "\n", sep="" )
      }
      VERBOSE( verbose, "nchoose[tot]:", exp(nperm.tot), "\n" )

      if ( log(nperm)>nperm.tot )
        warning( "more permutation than available" )
      
      for ( j in 1:length(idx) )
        perm[,idx[[j]]] <-
          t(apply(matrix(jar[[j]],nperm,length(jar[[j]]),byrow=T),1,sample))
    }
    return(perm)
  }
  # binary response variable
  #
  n1 <- sum(cls==lev[1])
  n2 <- sum(cls==lev[2])
  i.max <- which.max(c(n1,n2))
  max.lev <- lev[i.max]
  min.lev <- lev[-i.max]
  idx.max <- (1:m)[cls==max.lev]
  idx.min <- (1:m)[cls==min.lev]
  n.max <- c(n1,n2)[i.max]
  n.min <- c(n1,n2)[-i.max]

  if ( balanced )
  {
    if ( exhaustive ) {
      stop( "exhaustive not implemented for balanced permutations yet" )
    }
    perm[,idx.min] <- min.lev
    if ( length(lev)>2 )
      stop( "balanced permutation requires binary class" )
    if ( n1!=n2 && !equalized )
      stop( "balanced permutation requires same-size classes (or 'equalized=T')")
    k1 <- n.min%/%2
    k2 <- n.min%%2
        
    for ( i in 1:nperm )
    {
      k <- k1 + rbinom(1,1,.5)*k2
      idx.max.perm <- if ( n1!=n2 ) sample(idx.max,n.min) else idx.max
      perm[i,idx.max.perm] <- max.lev
      idx.max.perm <- sample( idx.max.perm, k )
      idx.min.perm <- sample( idx.min, 2*k-k )
      perm[i,idx.max.perm] <- min.lev
      perm[i,idx.min.perm] <- max.lev
    }
  }
  # IF NOT controlling for confounding variables
  #
  else if ( is.null(control) )
  {
    if ( exhaustive )
    {
      if ( choose(m,n.max)<=nperm )
      {
        nperm <- round(choose(m,n.max))
        VERBOSE( verbose, "Number of exhaustive permutations:", nperm, "\n" )
        nidx <- t(combn(m, n.max))
        perm  <- matrix( min.lev, nrow=nrow(nidx), ncol=m )
        for ( i in 1:nrow(nidx) ) {
          perm[i,nidx[i,]] <- max.lev
        }
        return( perm )
      }
      else {
        VERBOSE( T, "Number of exhaustive permutations is too large:",
                 choose(m,n.max),"(ignored)\n")
      }
    }
    # ELSE (if not exhaustive) ..
    #
    if ( log(nperm)>lchoose(n1+n2,n1)-log(2) )
      warning( paste("Number of possible permutations less than needed (", nperm, "): ",
                     choose(n1+n2,n1)/2, sep="") )
    for ( i in (1:nperm) )
    {
      idx <- if (equalized) 
               c(idx.min, sample(idx.max,n.min))
             else
               1:m
      perm[i,idx] <- sample(cls[idx],length(idx))
    }
  }
  # ELSE (controlling for confounding covariates) ..
  #
  else {
    lev0 <- lev[1]
    lev1 <- lev[2]
    control <- many2one(control)
    levs <- sort(unique(control))
    idx <- lapply( levs, function(z){which(control==z)} )
    nctl <- tabulate( control[cls==lev1], nbins=nlevels(control) )

    nperm.tot <- 0
    VERBOSE( verbose, "\n\tcontrolling for ", length(levs), "-level confounder\n", sep="" )

    # check if enough permutations available
    #
    for ( i in 1:length(nctl) ) {
      nperm.tot <- nperm.tot + (tmp <- lchoose(length(idx[[i]]),nctl[i]))
      VERBOSE( verbose, "\t  nchoose[",levels(control)[i],"]:   ", exp(tmp), "\n", sep="" )
    }
    VERBOSE( verbose, "\tnchoose[tot]:", exp(nperm.tot), "\n" )
    
    idx[nctl==0] <- NULL
    nctl <- nctl[nctl>0]
    if (nperm.tot<log(nperm)) {
      VERBOSE(verbose,"\tNumber of possible permutations less than required, reducing.\n" )
      nperm <- exp(nperm.tot)
      perm <- perm[1:nperm,]
    }    
    for ( i in 1:length(nctl) ) idx[[i]] <- c( nctl[i], idx[[i]] )
    
    for ( i in 1:nperm )
    {
      pidx <- unlist(sapply( idx, function(z){ sample( z[-1], z[1] ) } ))
      perm[i,] <- lev0
      perm[i,pidx] <- lev1
    }
  }
  return(perm)
}
