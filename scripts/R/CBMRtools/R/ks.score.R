# source( "~/dvlp/R/misc.R" )
# source( "~/dvlp/R/write.matrix.R" )
# source( "~/dvlp/R/levels.R" )
# source( "~/dvlp/R/permute.array.R" )
# source( "~/dvlp/R/read.res.R" )
# source( "~/dvlp/R/gene.cluster.R" )
# source( "~/dvlp/R/" )

# FUNCTION: KS GENESCORE
#
ks.genescore <- function
(
 n.x,               # length of ranked list
 y,                 # positions of geneset items in ranked list (basically, ranks)
 do.pval=T,         # compute asymptotic p-value
 alternative=c("two.sided","greater","less"),
 do.plot=F,         # draw the ES plot
 bare=F,            # return score & p-value only (a 2-tuple)
 weight=NULL,       # weights for weighted score (see Subramanian et al.) (usually, sort(score))
 weight.p=1,        # weights' exponent
 cls.lev=c(0,1),    # class labels to display
 absolute=F,        # takes max - min score rather than the maximum deviation from null
 plot.labels=FALSE, # hits' labels
 exact=NULL,
 ...                # additional plot arguments
 )
{
  # efficient version of ks.score (should give same results as ks.test, when weight=NULL)
  #
  alternative <- match.arg(alternative)
  DNAME <- paste( "1:", n.x, " and ", deparse(substitute(y)), sep="" )
  METHOD <- "Two-sample Kolmogorov-Smirnov test"
  n.y <- length(y)
  if ( n.y < 1 )  stop("Not enough y data")
  if ( any(y>n.x) ) stop( "y must be <= n.x: ", max(y) )
  if ( any(y<1) ) stop( "y must be positive: ", min(y) )
  if ( do.pval && !is.null(weight) ) warning("p-value meaningless w/ weighted score")
  if ( !is.null(weight) && length(weight)!=n.x ) stop("weights must be same length as ranked list: ", length(weight), " vs ", n.x)
  x.axis <- y.axis <- NULL

  # weighted GSEA score
  #
  if ( !is.null(weight) )
  {
    weight <- abs(weight[y])^weight.p

    Pmis <- rep(1, n.x); Pmis[y] <- 0; Pmis <- cumsum(Pmis); Pmis <- Pmis/(n.x-n.y)
    Phit <- rep(0, n.x); Phit[y] <- weight; Phit <- cumsum(Phit); Phit <- Phit/Phit[n.x]
    z <- Phit-Pmis

    score <- if (absolute) max(z)-min(z) else z[which.max(abs(z))]
    
    x.axis <- 1:n.x
    y.axis <- z
  }
  # KS score
  #
  else
  {
    y <- sort(y)
    n <- n.x * n.y/(n.x + n.y)
    hit <- 1/n.y
    mis <- 1/n.x

    # to compute score, only the y positions and their immediate preceding
    # ..positions are needed
    #
    Y <- sort(c(y-1,y)); Y <- Y[diff(Y)!=0]; y.match <- match(y,Y); if ( any(is.na(y.match)) ) browser()
    D <- rep( 0, length(Y) ); D[y.match] <- (1:n.y)
    zero <- which(D==0)[-1]; D[zero] <- D[zero-1]

    z <- D*hit - Y*mis
    
    score <- if (absolute) max(z)-min(z) else z[which.max(abs(z))]

    if (do.plot) {
      x.axis <- Y;
      y.axis <- z;
      if(Y[1]>0) {
        x.axis <- c(0,x.axis);
        y.axis <- c(0,y.axis);
      }
      if ( max(Y)<n.x ) {
        x.axis <- c(x.axis,n.x)
        y.axis <- c(y.axis,0)
      }
    }
  }
  if ( do.plot )
  {
    plot( x.axis, y.axis, type="l",
          xlab=paste("up-regulated for class ", cls.lev[2], " (KS>0) vs ",
                     "up-regulated for class ", cls.lev[1], " (KS<0)", sep="" ),
          ylab="gene hits",...)
    abline(h=0)
    abline(v=n.x/2,lty=3)
    axis(1,at=y,labels=plot.labels,tcl=0.25,las=2)
    i.max <- which.max(abs(y.axis))
    points( x.axis[i.max], y.axis[i.max], pch=20, col="red")
    text(x.axis[i.max]+n.x/20,y.axis[i.max],round(y.axis[i.max],2))
  }
  if ( !do.pval )
    return(score)

  # ELSE, compute asymptotic p-value
  #
  names(score) <- switch(alternative, two.sided="D", greater="D^+", less="D^-")
  PVAL <- ks.test(1:n.x,y=y,alternative=alternative,exact=exact)$p.value
  
  if ( bare ) {
    return( c(score=score, p.value=PVAL) )
  }
  RVAL <- list(statistic = score,
               p.value = PVAL, alternative = alternative, 
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}
# FUNCTION: GSET 2 LIST
#
gset2list <- function( gset, verbose=F )
{
  # transfer genesets from matrix format to list format
  #
  if ( ncol(gset)>1 ) {
    tmp <- NULL
    VERBOSE(verbose,
            "\tMapping", ncol(gset), "genesets ..")
    for ( i in 1:ncol(gset) ) {
      tmp <- c(tmp,
               list(unique(gset[gset[,i]!="NULL" & gset[,i]!="null" & gset[,i]!="",i])))
    }
    names(tmp) <- colnames(gset)
    gset <- tmp
    VERBOSE(verbose,
            " number of entries in each:\n", sep="")
    if (verbose) print( quantile(sapply(gset, length)))
    VERBOSE(verbose, "\n")
  }
  else {
    gset <- unique(gset[,1])
    VERBOSE( verbose, " (", length(gset), " entries).\n", sep="" )
  }
  gset
}
# FUNCTION: GSET 2 IDX
#
glist2idx <- function( glist, names, min.gset )
{
  # replace list of geneset names with list of geneset indices
  #
  if ( is.list(glist) )
  {
    glist.idx <- lapply( glist, function(z){ z.idx <- match(z,names);
                                             z.idx <- z.idx[!is.na(z.idx)]})
    tmp <- tmp.names <- NULL
    for ( i in 1:length(glist.idx) ) {
      if ( length(glist.idx[[i]])>=min.gset ) {
        tmp <- c( tmp, list(glist.idx[[i]]) )
        tmp.names <- c( tmp.names, names(glist.idx)[i] )
      }
    }
    names(tmp) <- tmp.names
    glist.idx <- tmp
    if ( (tmp <- length(glist)-length(glist.idx))>0 )
      cat("  Removed ", tmp, " genesets because too short (<",min.gset,")\n",sep="")
    if ( length(glist.idx)==0 )
      stop( "no geneset with the required minimum of genetags present in the dataset" )
  }
  else
  {
    glist.idx <- match( glist, names ); glist.idx <- glist.idx[!is.na(glist.idx)]
    n.idx <- length(glist.idx)
    if ( n.idx<min.gset ) {
      stop( paste("less tags than allowed (min=",min.gset,"): ", n.idx, sep="") )
    }
    glist.idx <- list(geneset=glist.idx)
  }
  glist.idx
}
# FUNCTION: KS PERM
#
ks.perm <- function
(
 x,                        # an (m-genes x n-samples) matrix
 cls=NULL,                 # class template
 tag=NULL,                 # id of gene to use as template (nearest neighbor mode)
 gset,                     # a vector or a matrix w/ one geneset per column, or a list of vectors
 min.gset=5,               # minimum length of accepted genesets
 robust=F,                 # use median instead of mean
 nperm=100,                # number of permutation iterations
 score=c("t.score","snr"), # differential score used when cls is specified
 paired=F,                 # paired differntial score?
 method=c("pearson","spearman","euclidean"),
                           # score used when tag is specified
 alternative=c("two.sided","greater","less"),
 weighted=F,               # use weighted KS score as in GSEA
 weight.p=1,               # weight's exponent
 smoother=0,               # smooth p-values
 control=NULL,
 cls.lev=c(0,1),
 seed=NULL, 
 do.pval=!weighted,
 balanced=F,
 exhaustive=F,
 verbose=F,
 do.plot=F,
 plot.name=NULL,
 plot.dev=pdf,
 do.abs=F,
 do.debug=F,
 ...
 )
{
  # compute ks-based enrichment wrt to class template (cls) or wrt
  # distance from gene (tag)
  #
  
  # some checks on inputs
  #
  if (is.null(dim(x))) stop( "x must be a 2D matrix" )
  if (is.null(cls) && is.null(tag) ) stop( "must specify either cls or tag" )
  if (!is.null(cls) && ncol(x)!=length(cls)) stop( "cls must be same length as ncol(x)" )
  #if (nperm<1 && !do.pval) stop( "must ask for either asymptotic or empirical p-value" )
  if (is.null(cls) && is.na(match(tag,rownames(x))) ) stop( "tag not found: ", tag )
  if (!is.null(dim(gset))) gset <- gset2list(gset, verbose=verbose)
  if (is.list(gset) && is.null(names(gset))) stop( "genesets must be named" )
  if (weighted && do.pval) stop("cannot compute asymptotic p-val w/ weighted score")
  if (!is.null(control) && is.null(levels(control)) ) stop( "control must have non-null levels" )
  
  ngset <- if (is.list(gset)) length(gset) else 1
  score  <- match.arg( score )
  method <- match.arg( method )
  alternative <- match.arg( alternative)

  # MANY-valued cls 
  #
  if ( !is.null(cls) && length(levels(cls))>2 )
  {
    levs <- levels(cls)
    KS <- NULL

    if ( !is.null(plot.name) ) {
      plot.dev( plot.name )
    }
    for ( lev in levs )
    {
      VERBOSE( verbose, "  ", lev, " vs. NOT ", lev, "\n", sep="")
      cls01 <- as.numeric( cls==match(lev,levs)-1 )
      levels(cls01) <- paste( c("not.",""), lev, sep="" )
      
      ks <- ks.perm(x,cls=cls01,gset=gset,min.gset=min.gset,robust=robust,nperm=nperm,
                    score=score,method=method,alternative=alternative,control=control,
                    weighted=weighted,weight.p=weight.p,cls.lev=levels(cls01),seed=seed,
                    do.pval=do.pval,exhaustive=exhaustive,verbose=verbose,
                    do.plot=do.plot || !is.null(plot.name),plot.name=NULL,plot.dev=plot.dev,
                    do.abs=do.abs,...)
      
      KS <- if (is.null(KS))
        ks
      else {
        cbind(KS, ks[match(rownames(KS),rownames(ks)),])
      }
    }
    if ( !is.null(plot.name) ) {
      dev.off()
      VERBOSE( verbose, "  plots saved to file '", plot.name, "'.\n", sep="" )
    }
    # yuck!!
    #
    if (nperm>0) {
      del <- mmatch( c("p1","fpr","GC.p","fdr.test"), colnames(KS) )
      del <- del[!is.na(del)]
      if (length(del)>0) KS <- KS[,-del]
      del <- which(colnames(KS)=="hits")
      KS <- KS[,-del[-length(del)]]
      del <- which(colnames(KS)=="size")
      KS <- KS[,-del[-length(del)]]
      colnames(KS)[which(colnames(KS)=="score")] <- paste( "score",levels(cls),sep=".")
    }
    else {
      if ( do.pval)
        colnames(KS)[which(colnames(KS)=="ks")] <- levels(cls)
      else
        colnames(KS) <- levels(cls)
    }
    return(KS)
  }
  # ELSE binary cls ...
  #
  dir <- if (alternative=="greater") 1 else -1

  SCORE <- if ( is.null(tag) )
    match.fun(score)
  else
    switch( method,
            euclidean = function(x,y){sqrt(sum((x-y)^2))},
            pearson =   function(x,y){1 - cor(x,y)},
            spearman =  function(x,y){1 - cor(rank(x),rank(y))} )

  ref.dist <- 1:nrow(x) # reference distribution
  gset.idx <- NULL
  ks.obs   <- ks.nul <- ks.pval <- results <- NULL
  gset.idx <- glist2idx( gset, names=rownames(x), min.gset=min.gset )
  ghit.len <- sapply( gset.idx, length )
  gset.len <- sapply( gset[match(names(gset.idx),names(gset))], length )
  ks.obs <- ks.pval <- rep( NA, length(gset.idx) );
  if (length(ks.obs)>1) names(ks.obs) <- names(gset.idx)

  # compute observed ranking scores
  #
  VERBOSE( verbose, "  computing observed rankings " )
  x.score <- if ( is.null(cls) ) {
    VERBOSE( verbose, "(based on NN) .. ")
    apply(x,1,SCORE,y=x[match(tag,rownames(x)),])
  }
  else {
    VERBOSE( verbose, "(based on class template) .. ")
    SCORE(x, cls=cls, robust=robust, paired=paired)
  }
  if (do.abs) {
    x.score <- abs(x.score)
  }
  x.rank <- rank(x.score,ties.method="first")
  VERBOSE( verbose, "done.\n" )
  VERBOSE( verbose, "  computing observed enrichment score(s) .." )
  if ( !is.null(plot.name) ) {
    plot.dev( plot.name )
  }
  for ( i in 1:length(gset.idx) )
  {
    tmp <- ks.genescore(nrow(x),y=x.rank[gset.idx[[i]]],do.pval=do.pval, bare=T, 
                        weight=if (weighted) sort(x.score), weight.p=weight.p,
                        do.plot=do.plot || !is.null(plot.name), cls.lev=cls.lev,
                        main=names(gset.idx)[i], ... )
    if ( do.pval ) {
      ks.obs[i] <- tmp[1]
      ks.pval[i] <- tmp[2]
    }
    else
      ks.obs[i] <- tmp
  }
  results <- { if (do.pval)
                 cbind(ks=ks.obs, asymptotic.p=ks.pval)
               else
                 matrix(ks.obs,length(ks.obs),1,dimnames=list(names(ks.obs),"ks"))
             }
  if ( !is.null(plot.name) ) {
    dev.off()
    VERBOSE( verbose, "(plots saved to file '", plot.name, "')\n", sep="" )
  }
  VERBOSE( verbose, " done.\n" )
  
  if (nperm<1)
    return(results)
  
  if ( !is.null(seed) ) set.seed(seed)

  # compute permuted ranking scores
  #  cls.perm has a permuted class labeling per row
  #
  VERBOSE( verbose, "  computing empirical p.value .. " )

  VERBOSE( verbose, "\n\tgenerating permuted labels .. " )

  cls.perm <- NULL
  if ( is.null(cls) )
  {
    y <- x[match(tag,rownames(x)),]
    cls.perm <- matrix( NA, nperm, length(y) )
    for ( i in 1:nperm ) {
      cls.perm[i,] <- sample( y, size=length(y) )
    }
  }
  else if (paired) {
    cls.perm <- permute.paired(cls,nperm=nperm,exhaustive=exhaustive,balanced=balanced)
  }
  else {
    cls.perm <- permute.binarray( cls, nperm=nperm, control=control, exhaustive=exhaustive,
                                  balanced=balanced, verbose=verbose)
    nperm <- nrow(cls.perm)
  }
  #VERBOSE( verbose, "done.\n" )

  x.prm <- matrix( 0, length(gset.idx), length(perm.2side.names) )
  colnames(x.prm) <- perm.2side.names
  
  # compute permuted ks scores
  #
  VERBOSE( verbose, "\tcomputing permuted ks scores .. " )
  percent <- pctstep <- max( .1, round(1/nperm,2) )
  for ( i in 1:nperm )
  {
    x.scores <- if (is.null(cls))
      apply(x,1,SCORE,y=cls.perm[i,])
    else if (paired)
      SCORE(x[,cls.perm[i,],drop=F],cls=cls,robust=robust,paired=T)
    else
      SCORE(x, cls=cls.perm[i,], robust=robust)
    if (do.abs) {
      x.scores <- abs(x.scores)
    }
    x.rnks <- rank(x.scores,ties.method="first")
    weight <- if (weighted) sort(x.scores)
    ks.nul <- sapply(gset.idx, function(z)
                     ks.genescore(nrow(x),y=x.rnks[z],weight=weight,weight.p=weight.p,do.pval=F) )

    #tmp <- perm.2side.online(ks.obs, ks.nul, dir=dir, do.abs=do.abs)
    tmp <- perm.2side.online(ks.obs, ks.nul, dir=dir)
    if ( !all(dim(tmp)==dim(x.prm)) ) stop("incompatible dimensions for summary")
    x.prm <- x.prm + tmp

    if ( verbose & i>=nperm*percent ) {
      VERBOSE( verbose, percent*100,"% ", sep="" )
      percent <- round(percent + pctstep,1)
    }
  }
  VERBOSE( verbose, "\n" )

  # complete computation of some p-values (see perm.2side.summary 
  # ..for the "offline" computation of these statistics
  #
  #x.out <- p2ss.add( x.prm=x.prm, x.obs=ks.obs, nperm=nperm, do.abs=do.abs )
  x.out <- p2ss.add( x.prm=x.prm, x.obs=ks.obs, nperm=nperm, smoother=smoother )

  bind.fun <- if (ngset>1) cbind else c
  names(ghit.len) <- names(gset.len) <- NULL
  if ( do.pval ) x.out <- bind.fun(x.out,
                                   asymptotic.p=signif(ks.pval),
                                   asymptotic.fdr=signif(pval2fdr(ks.pval)))
  x.out <- bind.fun( x.out, hits=ghit.len, size=gset.len )

  if ( do.debug )
    return( list(sum=x.out,prm=x.prm) )
  else
    return( x.out )
}
# FUNCTION: KS WRAPPER
#
ks.wrapper <- function( gp.filename, cls.filename, gset.filename, outstub=NULL, 
                        score=c("snr","t.score"), paired=F, weighted=F, weight.p=1,
                        min.gset=5, robust=F, nperm=0, tag.skip=0, do.pval=!weighted,
                        do.plot=!is.null(outstub), plot.dev=c("pdf","jpeg"),
                        verbose=F, header=T )
{
  # checks on inputs
  #
  if ( is.null(outstub) & (do.plot) )
    stop( "must specify oustub when asking for a plot" )

  plot.name <- if (do.plot) paste( outstub, ".", plot.dev, sep="" ) else NULL
  plot.dev <- match.arg( plot.dev )
  DEV <- match.fun( plot.dev )
  score <- match.arg(score)
  #score <- match.fun(score)
  file.test( gp.filename )
  if ( !is.null(cls.filename) ) file.test( cls.filename )
  file.test( gset.filename )

  # automatically load combinat from CRAN if necessary
  #
  #if (!require("combinat", quietly=TRUE)) {
  #  VERBOSE( verbose, "Installing package 'combinat'.\n" )
  #  install.packages("combinat")
  #}
  res <- read.data.gp( gp.filename )
  res <- res@signal
  cls  <- read.cls( cls.filename, verbose=verbose )
  VERBOSE( verbose, "  Reading genesets .. " )

  if ( length(grep(".Rout$",gset.filename,perl=T))>0 )
  {
    VERBOSE(verbose,"('.Rout' format) ")
    load(file=gset.filename)
    if (!exists("gset")) stop("no object named 'gset' was loaded")
  }
  else
  {
    is.gmt <- length(grep(".gmt",gset.filename))>0
    gset <- read.table( gset.filename, fill=T, header=!is.gmt, skip=tag.skip, sep="\t" )
    if ( is.gmt ) {
      VERBOSE( verbose, "(gmt format)" )
      gset <- gset[,-2] # eliminate the color column
      gset <- t(gset)
      colnames(gset) <- gset[1,]; gset <- gset[-1,,drop=F]
    }
    else {
      VERBOSE( verbose, "(gmx format)" )
      gset <- gset[-1,,drop=F] # eliminate the color row
    }
    if ( !header ) {
      colnames(gset) <- paste( "geneset.", 1:ncol(gset), sep="" )
    }
    gset <- gset2list( gset, verbose=verbose )  
  }
  VERBOSE( verbose, " done,", length(gset), "found.\n" )
  
  if ( ncol(res)!=length(cls) ) stop( "cls must be same length as ncol(res)" )

  results <- ks.perm( x=res, cls=cls, gset=gset, min.gset=min.gset, score=score,
                      paired=paired, weighted=weighted, weight.p=weight.p, robust=robust,
                      nperm=nperm, do.pval=do.pval, plot.name=plot.name, 
                      plot.dev=DEV, verbose=verbose, cls.lev=levels(cls) )
  
  rm.idx <- match( c("p1","GC.p","fpr","fdr.test"), colnames(results) )
  rm.idx <- rm.idx[!is.na(rm.idx)]
  results <- round(results[,-rm.idx,drop=F],4)

  if ( !is.null(outstub) ) {
    outname <- paste(outstub, ".txt", sep="")
    write.table( results, sep="\t", row.names=T, file=outname, quote=F )
    VERBOSE( verbose, "  Results saved to file '", outname, "'.\n\n", sep="" )
  }
  results
}
if (F) {  
CBMMLAB <- Sys.getenv('CBMMLAB')
  if (CBMMLAB=="") stop( "Use 'setenv CBMMLAB ..' to set CBMrepository's base directory" )
  RHOME <- paste(CBMMLAB,"R",sep="/")
  
  source( paste(RHOME, "source.ks.score.R", sep="/") )
  source( paste(RHOME, "hyper.enrichment.R", sep="/") )

  ## DAT is a [gene x sample] matrix
  ## rownames(DAT) must be the gene symbols

  DAT <- read.table(...)

  ## uplodad the geneset compendium of interest ..
  ## .. (and format it as a list, and turn all names to uppercase)
  ##
  GSETS <- select.enrichment.db("MSigDB_c2.cp",vers="v3.1")
  rownames(GSETS) <- toupper(rownames(GSETS))
  GSETS.list <- apply(GSETS,2,function(z) rownames(GSETS)[z==1])

  ## run GSEA
  ##
  GSEA <- ks.perm(DAT,cls=CLS,control=CTL,score="t.score",min.gset=5,gset=GSETS.list,
                  weighted=T,smoother=1,nperm=1000,verbose=T)
}
