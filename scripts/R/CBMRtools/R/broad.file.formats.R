# source("~/dvlp/R/read.res.new.R")

## simple classes to support Broad's '.res' and '.gct' format
##
## this is not well defined. Ideally, there would be overload of
## operators such as rownames, colnames, and indexing ([,]) so as to
## make it transparent. 

## use slotNames, getSlots, or getClass to see details about the class

###############################################################################################
#' gctdata
###############################################################################################
#'
#' Class \code{gctdata} represents a gene expression microarray data in the Broad '.gct' format
#'
#' @slot signal is the [m x n] gene-by-sample matrix of expression levels
#' @slot description is an m-sized character vector of genes' (or probesets') descriptions
#' 
#' @export
#' 
setClass("gctdata",
         representation=representation(
             signal="matrix",
             description="character"),
         validity=function(object) {
             if ( is.null(object@signal) )
                 stop( "is.null(signal)" )
             if ( is.null(object@description) )
                 stop( "is.null(description)" )
             if ( nrow(object@signal)!=length(object@description) )
                 stop('nrow(signal)!=length(description)')
             return(TRUE)
         })

## gctdata constructor

setMethod("initialize", "gctdata",
  function(.Object, signal=matrix(), description=character(1)) {
      .Object@signal <- signal
      .Object@description <- description
      #cat("*** Class: gctdata ***\n")
      #cat(">>> signal:      [", paste(dim(.Object@signal),collapse=' x '), "] matrix\n", sep="")
      #cat(">>> description: ", length(.Object@description), "-sized character vector\n", sep="")
      validObject(.Object)
      return(.Object)
  })

## GET
setGeneric("getSignal",function(object){standardGeneric("getSignal")})
#' @export
setMethod("getSignal","gctdata", function(object) { return(object@signal) })
#' @export
setGeneric("getDescription",function(object){standardGeneric("getDescription")})
setMethod("getDescription","gctdata", function(object) { return(object@description) })

## SET SIGNAL
setGeneric("setSignal<-",function(object,value){standardGeneric("setSignal<-")})
#' @export
setReplaceMethod(f="setSignal",
                 signature="gctdata",
                 definition=function(object,value) {
                     object@signal <- value
                     return(object)
                 })
## SET DESCRIPTION
setGeneric("setDescription<-",function(object,value){standardGeneric("setDescription<-")})
#' @export
setReplaceMethod(f="setDescription",
                 signature="gctdata",
                 definition=function(object,value) {
                     object@description <- value
                     return(object)
                 })

## (GET/SET) GENENAMES

setGeneric("genenames",function(object){standardGeneric("genenames")})
#' @export
setMethod(f="genenames",signature="gctdata",
          definition=function(object) { return ( rownames(object@signal) ) } )

setGeneric("genenames<-",function(object,value){standardGeneric("genenames<-")})
#' @export
setReplaceMethod(f="genenames",
                 signature="gctdata",
                 definition=function(object,value) {
                     rownames(object@signal) <- value
                     return ( object )
                 })

## (GET/SET) EXPTNAMES

setGeneric("exptnames",function(object){standardGeneric("exptnames")})
#' @export
setMethod(f="exptnames",signature="gctdata",
          definition=function(object) { return ( colnames(object@signal) ) } )

setGeneric("exptnames<-",function(object,value){standardGeneric("exptnames<-")})
#' @export
setReplaceMethod(f="exptnames",
                 signature="gctdata",
                 definition=function(object,value) {
                     colnames(object@signal) <- value
                     return ( object )
                 })

## OPERATOR '['

## READ GCT
##
read.gct <- function( file, force.read=FALSE, do.save=TRUE, binext=".RData",verbose=TRUE )
{
  # see if a binary object was cached
  #
  binfile <- gsub(".gct",binext,file)
  if ( !force.read && file.access(binfile)==0 ) {
    dat <- load.var(binfile)
    if (class(dat)!="gctdata") {
      stop("a gctdat object expected in '", binfile, "'", sep="")
    }
    return(dat)
  }
  VERBOSE( verbose, "\tReading signal .. " )
  x <- read.delim(file,sep="\t",fill=TRUE,skip=2, header=TRUE,check.names=FALSE,comment.char="")
  rownames(x) <- x[,1]
  desc <- as.character(x[,2])
  x <- as.matrix(x[,-c(1,2)])
  VERBOSE( verbose, "done, ", nrow(x), " genes x ", ncol(x), " experiments.\n", sep="" )

  dat <- new( "gctdata", signal=x, description=desc )
  if (do.save) {
    VERBOSE(verbose, "Saving binary object ..")
    save(dat,file=binfile)
    VERBOSE(verbose, " done.\n" )
  }
  return(dat)
}
## WRITE GCT
##
write.gct <- function( x, file, binext=".RData", do.save=TRUE, verbose=FALSE )
{
  if ( is.matrix(x) ) {
    x <- new("gctdata",
             signal=x,
             description=rownames(x))             
  }
  else if ( class(x)!="gctdata" )
    stop( "object does not belong to gctdata class" )
  
  cat( "#1.2\n", file=file )
  cat( paste(dim(x),collapse="\t"), sep="" , "\n", file=file, append=TRUE )
  cat( paste( c("Name\tDescription",exptnames(x)),collapse="\t"), sep="", "\n",
       file=file, append=TRUE )
  my.write.table(cbind(genenames(x),x@description,x@signal), row.names=FALSE, col.names=FALSE,
                 sep="\t", file=file, append=TRUE )

  if ( do.save ) {
    VERBOSE(verbose," (saving binary object ..")
    binfile <- gsub("\\.gct",binext,file)
    if ( binfile==file )
      warning("didn't save binary file (problems creating output name)")
    else {
      dat <- x
      save(dat,file=binfile)
    }
    VERBOSE(verbose, " done)")
  }
}
## DIM GCT
##
#' @export
dim.gctdata <- function( x )
{
  if ( class(x)!="gctdata" )
    stop( "object does not belong to gctdata class" )
  dim(getSignal(x))
}
## backward compatibility
##
set.exptnames <- function( x, enames ) {
    exptnames(x) <- enames
    return(x)
}
set.descriptors <- function( x, dnames ) {
    setDescprition(x) <- dnames
    return(x)
}
## SUBSET GCTDATA
##
subset.gctdata <- function( x, exptnames=NULL, genenames=NULL, ignore=FALSE )
{
  if (class(x)!="gctdata") stop( "'gctdata' object expected" )
  if ( is.null(exptnames) && is.null(genenames) )
    stop( "either exptnames or genenames must be provided" )

  g.idx <- 1:nrow(x)
  e.idx <- 1:ncol(x)
  
  # sample subsetting
  #
  if ( !is.null(exptnames) )
  {
    e.idx <- match.nona( exptnames, exptnames(x) )
    if ( !ignore && any(is.na(e.idx)) )
      stop( "subsetting on non-existing columns" )
    else
      e.idx <- e.idx[!is.na(e.idx)]
    if ( length(e.idx)==0 )
      stop( "empty subset to select" )
  }
  # gene subsetting
  #
  if ( !is.null(genenames) )
  {
    g.idx <- match( genenames, genenames(x) )
    if ( !ignore && any(is.na(g.idx)) )
      stop( "subsetting on non-existing rows" )
    else
      g.idx <- g.idx[!is.na(g.idx)]
    if ( length(g.idx)==0 )
      stop( "empty subset to select" )
  }
  return( new("gctdata",
              signal=getSignal(x)[g.idx,e.idx,drop=FALSE],
              description=getDescription(x)[g.idx]) )
}
###############################################################################################
#' resdata
###############################################################################################
#'
#' Class \code{resdata} represents a gene expression microarray data in the Broad '.res' format
#'
#' @slot signal is the [m x n] gene-by-sample matrix of expression levels
#' @slot calls is the [m x n] gene-by-sample matrix of P/M/A calls
#' @slot description is an m-sized character vector of genes' (or probesets') descriptions
#' @slot scale is the scaling factor used with each sample (if any)
#' 
#' @export
#' 
setClass("resdata",
       contains="gctdata",
       representation(calls="matrix",
                      scale="character"),
       validity = function(object) {
           if ( is.null(object@calls) )
               stop( "is.null(calls)" )
           if ( ncol(object@signal)!=ncol(object@calls) )
               stop( "ncol(signal)!=ncol(calls)" )
           if ( nrow(object@signal)!=nrow(object@calls) )
               stop( "nrow(signal)!=nrow(calls)" )
           if ( any(colnames(object@signal)!=colnames(object@calls)) )
               stop( "colnames(signal)!=colnames(calls)" )
           if ( any(rownames(object@signal)!=rownames(object@calls)) )
               stop( "rownames(signal)!=rownames(calls)" )
           if ( !is.null(object@scale) && (length(object@scale)!=ncol(object@signal)) )
               stop( "length(scale)!=ncol(signal)" )
           return(TRUE)
       })

## resdata CONSTRUCTOR

setMethod("initialize", "resdata",
  function(.Object,signal=matrix(),calls=matrix(),scale=NULL,description=character(1)) {
    .Object@calls=calls
    .Object@scale=scale
    .Object <- callNextMethod(.Object,signal=signal,description=description)
    #cat("*** Class: resdata ***\n")
    #cat(">>> signal:      [", paste(dim(.Object@signal),collapse=' x '), "] matrix\n", sep="")
    #cat(">>> calls:       [", paste(dim(.Object@calls),collapse=' x '), "] matrix\n", sep="")
    #cat(">>> scale:       ", length(.Object@scale), "-sized character vector\n", sep="")
    #cat(">>> description: ", length(.Object@description), "-sized character vector\n", sep="")
    validObject(.Object)
    return(.Object)
  })

## READ RES

read.res <- function
(
 file,
 nrow=TRUE,
 sep="\t",
 ignore.scale=FALSE,
 description="Description",
 accession="Accession",
 force.read=FALSE,
 do.save=TRUE,
 binext=".RData",
 verbose=TRUE
 )
{
  # see if a binary object was cached
  #
  binfile <- gsub("\\.res",binext,file)
  if ( !force.read && file.access(binfile)==0 ) {
    VERBOSE( verbose, "Loadin bin file '", binfile, "' ..", sep="")
    dat <- load.var(binfile)
    if (class(dat)!="resdata") {
      stop("a resdat object expected in '", binfile, "'", sep="")
    }
    VERBOSE(verbose, "done, ",nrow(dat), " genes and ", ncol(dat), " experiments.\n", sep="")
    return(dat)
  }
  VERBOSE( verbose, "Reading file '", file, "' .. ", sep="" )

  calls <- desc <- NULL
  scale <- ""
  
  # Reading header line
  #
  chip.names <-
    sapply(read.delim(file,nrows=1,sep=sep,header=FALSE,colClasses="character",fill=TRUE),as.character)

  if (chip.names[length(chip.names)]=="") chip.names <- chip.names[-length(chip.names)]
  if (chip.names[1]!=description)
    stop("Unrecognized description label: ",chip.names[1]," (", description," expected)",sep="")
  if (chip.names[2]!=accession)
    stop("Unrecognized accession label: ",chip.names[1]," (", accession," expected)",sep="")
  
  # Reading data (signal + calls)
  #
  dat <- read.delim(file,skip=3,sep=sep,header=FALSE,quote="",check.names=FALSE,row.names=2)
  if ( length(chip.names)!=ncol(dat) )
    stop("header and dat have different length: ", length(chip.names), "!=", ncol(dat),sep="")

  # Formatting chip names
  #
  chip.names <- chip.names[seq(3,ncol(dat),2)]
  if ( any(na.idx <- is.na(chip.names)) )
    stop( "Chips' names are missing: ", paste(which(na.idx),collapse=","))

  # Reading scale factors
  #
  if ( !ignore.scale )
  {
    scale <- drop(as.matrix(read.delim(file,nrows=1,sep=sep,skip=1,header=FALSE)))[seq(3,ncol(dat),2)]
    names(scale) <- chip.names
    if ( any(na.idx <- is.na(scale)))
      stop( "Scale factors are missing: ", paste(which(na.idx),collapse=","))
  }
  # Reading number of genes
  #
  ngene <- read.delim(file,nrows=1,sep=sep,skip=2,header=FALSE)[1]
  if ( ngene!=nrow(dat) )
    stop("Number of rows different from expected: ", nrow(dat), " (",ngene,") expected", sep="")

  desc <- drop(as.matrix(dat[,1]))
  na.idx <- is.na(desc)
  desc[na.idx] <- "No information for gene"
  calls <- as.matrix(dat[,seq(3,ncol(dat),2)])
  dat <- as.matrix(dat[,seq(2,ncol(dat),2)])
  colnames(dat) <- colnames(calls) <- chip.names
  
  VERBOSE(verbose, "done, ",nrow(dat), " genes and ", ncol(dat), " experiments.\n", sep="")

  dat <- new("resdata",
             signal=dat,
             calls=calls,
             scale=scale,
             description=desc)
  if (do.save) {
    VERBOSE(verbose, "Saving binary object ..")
    save(dat,file=binfile)
    VERBOSE(verbose, " done.\n" )
  }
  dat
}
## DIM RESDATA
##
#' @export
dim.resdata <- function( x )
{
  if ( class(x)!="resdata" )
    stop( "object does not belong to resdata class" )
  dim(getSignal(x))
}
## MERGE RESDATA
##
merge.resdata <- function( X, Y )
{
  if (class(X)!="resdata") stop("resdata expected for X")
  if (class(Y)!="resdata") stop("resdata expected for Y")
  if ( length(expts <- intersect(exptnames(X),exptnames(Y)))>0 )
    stop("non unique experiment names: ", paste( expts, collapse=", " ), sep="")
  genes <- match(genenames(X),genenames(Y))
  if (any(is.na(genes)))   stop( "X and Y don't have the same genes" )

  new("resdata",
      signal=cbind(getSignal(X),getSignal(Y)[genes,]),
      calls =cbind(getCalls(X),getCalls(Y)[genes,]),
      description=getDescription(X),
      scale=c(getScale(X),getScale(Y)) )  
}
## SUBSET RESDATA
##
subset.resdata <- function( x, exptnames=NULL, genenames=NULL, ignore=FALSE )
{
  if (class(x)!="resdata") stop( "'resdata' object expected" )
  if ( is.null(exptnames) && is.null(genenames) )
    stop( "either exptnames or genenames must be provided" )

  g.idx <- 1:nrow(x)
  e.idx <- 1:ncol(x)
  
  # sample subsetting
  #
  if ( !is.null(exptnames) )
  {
    e.idx <- match( exptnames, colnames(getSignal(x)) )
    if ( !ignore && any(is.na(e.idx)) )
      stop( "subsetting on non-existing columns" )
    else
      e.idx <- e.idx[!is.na(e.idx)]
    if ( length(e.idx)==0 )
      stop( "empty subset to select" )
  }
  # gene subsetting
  #
  if ( !is.null(genenames) )
  {
    g.idx <- match( genenames, rownames(getSignal(x)) )
    if ( !ignore && any(is.na(g.idx)) )
      stop( "subsetting on non-existing rows" )
    else
      g.idx <- g.idx[!is.na(g.idx)]
    if ( length(g.idx)==0 )
      stop( "empty subset to select" )
  }
  return( new("resdata",
              signal=getSignal(x)[g.idx,e.idx,drop=FALSE],
              calls=getCalls(x)[g.idx,e.idx,drop=FALSE],
              description=getDescription(x)[g.idx],
              scale=if (length(getScale(x))>1) getScale(x)[e.idx] else "") )
}
## RENAME RESDATA
##
rename.data <- function( x, expt.names=NULL, gene.names=NULL, descriptors=NULL, new.names )
{
  ## EXPERIMENT name substitution
  ##
  if ( !is.null(expt.names) ) {
    if (length(expt.names)!=length(new.names))
      stop("exptnames and new.names must be same length")
    if ( any(is.na(match.idx <- match(expt.names,exptnames(x)))) )
      stop( "missing exptnames" )

    colnames(x@signal)[match.idx] <- new.names
    if (class(x)=="resdata") {
      colnames(x@calls)[match.idx] <- new.names
      names(x@scale)[match.idx] <- new.names
    }
    #if (class(x)=="l1000data") {
    #  colnames(x@meta)[match.idx] <- new.names
    #}
  }
  ## GENE name substitution
  ##
  else if ( !is.null(gene.names) )
  {
    if (length(gene.names)!=length(new.names))
      stop("exptnames and new.names must be same length")
    if ( any(is.na(match.idx <- match(gene.names,genenames(x)))) )
      stop( "missing genenames" )
    rownames(x@signal)[match.idx] <- new.names
    if (class(x)=="resdata") {
      rownames(x@calls)[match.idx] <- new.names
    }
  }
  ## DESCRIPTION name substitution
  ##
  else if ( !is.null(descriptors) )
  {
    if (length(descriptors)!=length(new.names))
      stop("descriptors and new.names must be same length")
    if ( any(is.na(match.idx <- match(descriptors,descriptors(x)))) )
      stop( "missing descriptors" )
    x@description[match.idx] <- new.names
  }
  else
    stop( "either expt.names or gene.names or descriptors must be non-null")
  x
}
## SUBSET DATA
##
subset.data <- function( x, exptnames=NULL, genenames=NULL, ignore=FALSE )
{
  if (class(x)=="gctdata") {
    return( subset.gctdata(x, exptnames=exptnames, genenames=genenames, ignore=ignore) )
  }
  else if (class(x)=="resdata") {
    return( subset.resdata(x, exptnames=exptnames, genenames=genenames, ignore=ignore) )
  }
  #else if ( class(x)=="l1000data" ) {
  #  return( subset.l1000data(x, exptnames=exptnames, genenames=genenames, ignore=ignore) )
  #}
  else stop( "unrecognized data format: ", class(x) )
}
## READ CLS
##
read.cls <- function( file, rowskip=2, do.lbls=(rowskip==2), gc.lbl=FALSE, verbose=TRUE )
{
  # INPUT:
  #   - file     filename
  #   - rowskip  rows to skip before actual class labels
  #   - do.lbls  whether the 2nd row contains the class names (row starting with '#')
  #   - gc.lbl   type of mapping from class names to class labels (default is first ..
  #              ..class name associated to lowest class label, etc.)
  #
  stats <- as.vector(as.matrix(read.table(file,nrows=1,header=FALSE)))
  if ( length(stats)!=3 )
    stop( "1st row of '.cls' file must have three (3) entries" )
  cls <- as.vector(as.matrix(read.table(file, skip=rowskip, header=FALSE)))
  if ( length(cls)!=stats[1] )
    stop( paste("Wrong number of labels: ", length(cls), " (", stats[1], " expected)", sep="") )
  if ( my.nlevels(cls)!=stats[2] )
    stop( paste("Wrong number of label IDs: ", my.nlevels(cls), " (", stats[2], " expected)", sep="") )
  if (do.lbls)
  {
    lbls <- sapply(read.table(file,skip=1,nrows=1,comment.char=""),as.character)[-1]

    if ( length(lbls)!=stats[2] )
      stop( "label IDs must be as many as labels' number in '.cls' file" )

    # label IDs always sorted from the one associated to lowest label to highest label
    #
    levels(cls) <- if ( gc.lbl )
      lbls[order(unique(cls))]
    else
      lbls

    VERBOSE( verbose, "class labels: ", paste(levels(cls),collapse=", "), "\n" )
  }
  return( cls )    
}
## WRITE CLS
##
write.cls <- function( cls, filen )
{
  cat( length(cls), length(levels(cls)), "1\n", file=filen )
  cat( "#", paste( levels(cls) ), sep=" ", file=filen, append=TRUE )
  cat( "\n", file=filen, append=TRUE )
  cat( paste( cls ), sep=" ", file=filen, append=TRUE )
  cat( "\n", file=filen, append=TRUE )
}
## SUBSET CLS
##
subset.cls <- function( cls, sub, do.renum=FALSE )
{
  levs <- if (is.null(levels(cls))) sort(unique(cls)) else levels(cls)
  cls1 <- cls[sub]
  cls1.lev <- levs[match(sort(unique(cls1)),sort(unique(cls)))]
  if (do.renum) cls1 <- match(cls1,unique(cls1))-1
  levels(cls1) <- cls1.lev
  cls1
}
## TAB 2 CLS
##
tab2cls <- function(tab,col.name=NULL,cnames=NULL,levs=NULL,na.rm=TRUE)
{
  if ( is.null(col.name) && !is.null(dim(tab)) )
    stop( "must provide col.name with 2D tab" )
  if ( !is.null(dim(tab)) && !(col.name %in% colnames(tab)) )
    stop( "non-existing column: ", col.name )
  if ( !is.null(dim(tab)) )
    tab <- tab[,col.name]
  if ( !is.null(cnames) && length(cnames)!=length(tab) )
    stop( "cnames must be same length as tab" )
  if ( !is.null(levs) && length(levs)!=length(unique(tab)) )
    stop( "levs must be as many as class values" )
  
  cls <- if (na.rm) match(tab,unique(tab[!is.na(tab)]))-1 else match(tab,unique(tab))-1
  levels(cls) <- if (is.null(levs)) {if (na.rm) unique(tab[!is.na(tab)]) else unique(tab)} else levs
  if ( !is.null(cnames) )
    names(cls) <- cnames

  cls
}
## MANY 2 ONE CLS
##
many2one <- function( cls ) many2one.cls(cls) # just for backward compatibility
many2one.cls <- function( cls )
{
  if ( is.null(dim(cls)) || ncol(cls)==1 ) {
    cls2 <- as.numeric(match(drop(cls),sort(unique(drop(cls)))))
    levels(cls2) <- levels(cls)
    return(cls2)
  }
  # ELSE ..
  #
  levs <- sort(apply(unique(cls),1,paste,collapse="."))
  cls.new <- as.numeric( match( apply(cls,1,paste,collapse="."), levs) )
  levels(cls.new) <- levs
  cls.new
}
## CREATE CLS
##
create.cls <- function( X, levs=NULL, nms=NULL, do.sort=TRUE )
{
  if ( !is.null(nms) && !is.null(names(X)) )
    stop( "names(X) is non-null, can't assign new names" )
  if ( !is.null(nms) )
    names(X) <- nms
  
  if (do.sort) {
    ord <- order(X)
    nms <- names(X)
    X <- X[ord]
    if (!is.null(nms))
      names(X) <- nms[ord]
  }
  levels(X) <- if (is.null(levs)) unique(sort(X)) else levs
  X
}
## SUBSET ALL
##
subset.all <- function( dat, cls, subfile=NULL, subset=NULL, verbose=FALSE )
{
  if ( is.null(subset) && is.null(subfile) )
    stop( "subset or subfile must be specified" )

  VERBOSE( verbose, "Extracting sample subset " )
  if ( !is.null(subfile) ) {
    VERBOSE( verbose, "from '", subfile, "' .. ", sep="" )
    subset <- read.table( subfile, sep="\t", header=FALSE )[,1]
  }
  s.idx <- match( subset, exptnames(dat) )
  if ( any(is.na(s.idx)) ) stop( "subset contains unmatched samples" )
  dat <- subset.data( dat, exptnames=subset )
  cls <- subset.cls( cls, s.idx )
  VERBOSE( verbose, "done, [", paste(dim(dat),collapse=" x "), "] data entries.\n", sep="" )
  list( data=dat, cls=cls )
}

###############################################################################################
###############################################################################################
if (FALSE)
{
    CBMGIT <- Sys.getenv('CBMGIT')
    if (CBMGIT=="") stop( "Use 'setenv CBMGIT ..' to set CBMgithub's base directory" )
    source( paste(CBMGIT, "scripts/R/CBMRtools/R/misc.R", sep="/") )
    source( paste(CBMGIT, "scripts/R/broad.file.formats.dvlp.R", sep="/") )

    setwd('~/Research/Projects/oralcancer/taz_yap_dpagt1/processed_data/')

    GCT <- read.gct('DUMMY.gct',force.read=TRUE)
    
    print(genenames(GCT))
    genenames(GCT) <- paste("G",1:nrow(GCT),sep="")
    print(genenames(GCT))

    print(exptnames(GCT)[1:10])
    exptnames(GCT) <- paste("E",1:ncol(GCT),sep="")
    print(exptnames(GCT)[1:10])
    

    RES <- read.res('DUMMY.res',force.read=TRUE)
}

