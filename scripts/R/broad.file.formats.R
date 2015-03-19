# source("~/dvlp/R/read.res.new.R")

## simple classes to support Broad's '.res' and '.gct' format
##
## this is not well defined. Ideally, there would be overload of
## operators such as rownames, colnames, and indexing ([,]) so as to
## make it transparent. 

setClass("resdata",
         representation(signal="matrix",
                        calls="matrix",
                        description="character",
                        scale="character"))
         
setClass("gctdata",
         representation(signal="matrix",
                        description="character"))

setClass("l1000data",
         representation(signal="matrix",
                        meta="matrix",
                        description="character"),
         validity = function(object) {
           if ( ncol(object@signal)!=ncol(object@meta) )
             stop( "ncol(object@signal)!=ncol(object@meta)" )
           if ( any(colnames(object@signal)!=colnames(object@meta)) )
             stop( "colnames(object@signal)!=colnames(object@meta)" )
         }
         )


# to get the names of the class slots, use
#
# names(attributes(tmp))
#
# where tmp is an instance of the class

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
  chip.names <- sapply(read.delim(file,nrows=1,sep=sep,header=FALSE,colClasses="character",fill=TRUE),as.character)

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
write.res <- function( x, file=NULL, binext=".RData", verbose=TRUE, do.save=TRUE )
{
  if ( class(x)!="resdata" )
    stop( "object does not belong to resdata class" )

  res <- matrix( NA, nrow(x@signal), ncol(x@signal)*2 )
  idx <- seq( 1, ncol(res)-1, 2 )
  res[,idx]  <- x@signal
  res[,-idx] <- x@calls
  
  cat( "Description\tAccession\t",
      paste( colnames(x@signal), collapse="\t\t", sep="" ),
      "\t\n", sep="", file=file )
  cat("\t\t",
      paste( x@scale, collapse="\t\t", sep="" ),
      "\t\n", sep="", file=file, append=TRUE )
  cat(dim(x)[1], "\n", sep="", file=file, append=TRUE )

  my.write.table(cbind(x@description, rownames(x@signal), res ),
                 sep="\t", col.names=FALSE, file=file, append=TRUE )
  if ( do.save ) {
    VERBOSE(verbose," (saving binary object ..")
    binfile <- gsub("\\.res",binext,file)
    if ( binfile==file )
      warning("didn't save binary file (problems creating output name)")
    else {
      dat <- x
      save(dat,file=binfile)
    }
    VERBOSE(verbose, " done)")
  }
}
dim.resdata <- function( x )
{
  if ( class(x)!="resdata" )
    stop( "object does not belong to resdata class" )
  dim(x@signal)
}
genenames <- function( x )
{
  if ( !(class(x) %in% c("resdata","gctdata","l1000data")) )
    stop( "'resdata|gctdata' object expected" )
  rownames(x@signal)
}
exptnames <- function( x )
{
  if ( !(class(x) %in% c("resdata","gctdata","l1000data")) )
    stop( "'resdata|gctdata|l1000data' object expected" )
  colnames(x@signal)
}
signal <- function( x )
{
  x@signal
}
descriptors <- function( x )
{
  if ( !(class(x) %in% c("resdata","gctdata","l1000data")) )
    stop( "'resdata|gctdata' object expected" )
  x@description
}
set.genenames <- function( x, gnames )
{
  if ( length(gnames)!=nrow(x) )
    stop( "number of gene names doesn't match number of data rows" )
  if ( length(gnames)!=length(unique(gnames)) )
    stop( "gene names must be unique" )
  rownames(x@signal) <- gnames
  if ( class(x)=="resdata" )
    rownames(x@calls) <- gnames
  x
}
set.exptnames <- function( x, enames )
{
  if ( length(enames)!=ncol(x) )
    stop( "number of experiment names doesn't match number of data columns" )
  if ( length(enames)!=length(unique(enames)) )
    stop( "experiment names must be unique" )
  colnames(x@signal) <- enames
  if ( class(x)=="resdata" )
    rownames(x@calls) <- enames
  x
}
set.descriptors <- function( x, dnames )
{
  if ( length(dnames)!=nrow(x) )
    stop( "number of descriptors doesn't match number of data rows" )
  x@description <- dnames
  x
}
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
    e.idx <- match( exptnames, colnames(x@signal) )
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
    g.idx <- match( genenames, rownames(x@signal) )
    if ( !ignore && any(is.na(g.idx)) )
      stop( "subsetting on non-existing rows" )
    else
      g.idx <- g.idx[!is.na(g.idx)]
    if ( length(g.idx)==0 )
      stop( "empty subset to select" )
  }
  return( new("resdata",
              signal=x@signal[g.idx,e.idx,drop=FALSE],
              calls=x@calls[g.idx,e.idx,drop=FALSE],
              description=x@description[g.idx],
              scale=if (length(x@scale)>1) x@scale[e.idx] else "") )
}
merge.resdata <- function( X, Y )
{
  if (class(X)!="resdata") stop("resdata expected for X")
  if (class(Y)!="resdata") stop("resdata expected for Y")
  if ( length(expts <- intersect(exptnames(X),exptnames(Y)))>0 )
    stop("non unique experiment names: ", paste( expts, collapse=", " ), sep="")
  genes <- match(genenames(X),genenames(Y))
  if (any(is.na(genes)))   stop( "X and Y don't have the same genes" )

  new("resdata",
      signal=cbind(X@signal,Y@signal[genes,]),
      calls =cbind(X@calls,Y@calls[genes,]),
      description=X@description,
      scale=c(X@scale,Y@scale) )  
}
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
dim.gctdata <- function( x )
{
  if ( class(x)!="gctdata" )
    stop( "object does not belong to gctdata class" )
  dim(x@signal)
}
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
    e.idx <- match( exptnames, colnames(x@signal) )
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
    g.idx <- match( genenames, rownames(x@signal) )
    if ( !ignore && any(is.na(g.idx)) )
      stop( "subsetting on non-existing rows" )
    else
      g.idx <- g.idx[!is.na(g.idx)]
    if ( length(g.idx)==0 )
      stop( "empty subset to select" )
  }
  return( new("gctdata",
              signal=x@signal[g.idx,e.idx,drop=FALSE],
              description=x@description[g.idx]) )
}
read.l1000 <- function( fname, force.read=FALSE, do.save=TRUE, binext=".RData",verbose=TRUE )
{
  gene.headings <- c("pr_gene_symbol","pr_gene_title")
  
  ## see if a binary object was cached
  ##
  binfile <- gsub("\\.gct",binext,fname)
  if ( !force.read && file.access(binfile)==0 ) {
    dat <- load.var(binfile)
    if (class(dat)!="l1000data") {
      stop("a l1000data object expected in '", binfile, "'", sep="")
    }
    return(dat)
  }
  heading <- as.vector(as.matrix(read.tab.delim(fname,header=FALSE,skip=1,nrows=1)))

  VERBOSE( verbose, "\tReading signal .. " )
  x <- read.delim(fname,sep="\t",fill=TRUE,skip=2, row.names=1,header=TRUE,check.names=FALSE,comment.char="",stringsAsFactors=FALSE)
  if ( nrow(x)!=heading[1]+heading[4] )
    stop( "unexpected # of rows:", heading[1]+heading[4], "expected" )
  if ( ncol(x)!=heading[2]+heading[3] )
    stop( "unexpected # of columns:", heading[2]+heading[3], "expected" )
  if ( any(is.na(desc.idx <- match(gene.headings,colnames(x)))) )
    stop( "missing columns: ", gene.headings[is.na(desc.idx)] )
  
  meta <- x[1:heading[4],-(1:heading[3])]
  desc <- apply(x[,desc.idx],1,paste,collapse=":")[-(1:heading[4])]
  x <- as.matrix(x[-(1:heading[4]),-(1:heading[3])])
  rnames <- rownames(x)
  x <- apply(x,2,as.numeric)
  rownames(x) <- rnames
  dat <- new( "l1000data", signal=as.matrix(x), description=desc, meta=as.matrix(meta) )

  VERBOSE( verbose, "done, ", nrow(dat), " genes x ", ncol(dat), " experiments.\n", sep="" )

  if (do.save) {
    VERBOSE(verbose, "Saving binary object ..")
    save(dat,file=binfile)
    VERBOSE(verbose, " done.\n" )
  }
  return(dat)
}
dim.l1000data <- function( x )
{
  if ( class(x)!="l1000data" )
    stop( "object does not belong to l1000data class" )
  dim(x@signal)
}
subset.l1000data <- function( x, exptnames=NULL, genenames=NULL, ignore=FALSE )
{
  if (class(x)!="l1000data") stop( "'l1000data' object expected" )
  if ( is.null(exptnames) && is.null(genenames) )
    stop( "either exptnames or genenames must be provided" )

  g.idx <- 1:nrow(x)
  e.idx <- 1:ncol(x)
  
  # sample subsetting
  #
  if ( !is.null(exptnames) )
  {
    e.idx <- match( exptnames, colnames(x@signal) )
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
    g.idx <- match( genenames, rownames(x@signal) )
    if ( !ignore && any(is.na(g.idx)) )
      stop( "subsetting on non-existing rows" )
    else
      g.idx <- g.idx[!is.na(g.idx)]
    if ( length(g.idx)==0 )
      stop( "empty subset to select" )
  }
  return( new("l1000data",
              signal=x@signal[g.idx,e.idx,drop=FALSE],
              meta=x@meta[,e.idx,drop=FALSE],
              description=x@description[g.idx]) )
}
l1000.from.gctx <- function( GCTX, cnames=c("cell_id","pert_time","pert_id") )
{
  ## checks on inputs
  ##
  if ( !all(names(GCTX) %in% c("matrix","row.annotations","column.annotations")) )
    stop( "GCTX expected to have fields 'matrix','row.annotations','column.annotations'" )
  if ( !all(cnames %in% colnames(GCTX$column.annotations)) )
    stop( "cnames missing from GCTX$column.annotations' column names" )
  if ( !("id" %in% colnames(GCTX$row.annotations)) )
    stop( "'id' not found among row.annotations' column names" )
  if ( ncol(GCTX$matrix)!=nrow(GCTX$column.annotations) )
    stop( "ncol(GCTX$matrix)!=nrow(GCTX$column.annotations)" )
  ## end checks
  
  rownames(GCTX$matrix) <- GCTX$row.annotations[,'id']
  META <- t(as.matrix(GCTX$column.annotations))
  colnames(META) <- colnames(GCTX$matrix) <- apply(GCTX$column.annotations[,cnames],1,paste,collapse="_")
  
  DAT <- new("l1000data",
             signal=GCTX$matrix,
             meta=META,
             description=apply(GCTX$row.annotations[,c('pr_gene_symbol','pr_gene_title')],1,paste,collapse=" : ")) 
}
merge.data <- function( D1, D2 )
{
  if ( class(D2)!=class(D1) )
    stop( "D2's class must be same as D1's class:",class(D2) )
  if ( any(genenames(D1)!=genenames(D2)) )
    stop( "genenames(D1)!=genenames(D2)" )
  
  if ( class(D1)=="l1000data" )
  {
    if ( any(rownames(D1@meta)!=rownames(D2@meta)) )
      stop( "rownames(D1@meta)!=rownames(D2@meta)" )
    return(new("l1000data",
               signal=cbind(D1@signal,D2@signal),
               description=D1@description,
               meta=cbind(D1@meta,D2@meta)))
  }
  else if ( class(D1)=="gctdata" )
  {
    return(new("gctdata",
               signal=cbind(D1@signal,D2@signal),
               description=D1@description))
  }
  else if ( class(D1)=="resdata" )
  {
    return(new("resdata",
               signal=cbind(D1@signal,D2@signal),
               calls=cbind(D1@calls,D2@calls),
               scale=c(D1@scale,D2@scale),
               description=D1@description))
  }
  else {
    stop( "merge.data not implemented for class:",class(D1) )
  }
}
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
    if (class(x)=="l1000data") {
      colnames(x@meta)[match.idx] <- new.names
    }
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
subset.data <- function( x, exptnames=NULL, genenames=NULL, ignore=FALSE )
{
  if (class(x)=="gctdata") {
    return( subset.gctdata(x, exptnames=exptnames, genenames=genenames, ignore=ignore) )
  }
  else if (class(x)=="resdata") {
    return( subset.resdata(x, exptnames=exptnames, genenames=genenames, ignore=ignore) )
  }
  else if ( class(x)=="l1000data" ) {
    return( subset.l1000data(x, exptnames=exptnames, genenames=genenames, ignore=ignore) )
  }
  else stop( "unrecognized data format" )
}
read.data.gp <- function( data.filename, verbose=TRUE, ... )
{
  #file.test( data.filename )
  ext <- file.ext( data.filename )

  VERBOSE( verbose, "File extension: '.", ext, "'.\n", sep="" )

  res <- switch( ext,
                 res=read.res( data.filename, verbose=verbose, ... ),
                 gct=read.gct( data.filename, verbose=verbose, ... ),
                 odf=stop( "odf reader not implemented yet" ),
                 stop("unrecognized file extension: ", ext) )
}
write.data.gp <- function( data, filename, do.save=FALSE, verbose=TRUE )
{
  #file.test( filename, mode=2 )
  ext <- file.ext(filename)

  VERBOSE(verbose, "Writing data to file '", filename, "' ..", sep="");
  if ( class(data)=="resdata" ) {
    if ( ext!="res" ) warning( "improper extension for res file: ", ext )
    write.res( data, filename, do.save=do.save, verbose=verbose )
  }
  else if ( class(data)=="gctdata" ) {
    if ( ext!="gct" ) warning( "improper extension for gct file: ", ext )
    write.gct( data, filename, do.save=do.save, verbose=verbose )
  }
  else {
    stop( "unsupported data format: ", class(data) )
  }
  VERBOSE(verbose, "done.\n" )
}
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
write.cls <- function( cls, filen )
{
  cat( length(cls), length(levels(cls)), "1\n", file=filen )
  cat( "#", paste( levels(cls) ), sep=" ", file=filen, append=TRUE )
  cat( "\n", file=filen, append=TRUE )
  cat( paste( cls ), sep=" ", file=filen, append=TRUE )
  cat( "\n", file=filen, append=TRUE )
}
subset.cls <- function( cls, sub, do.renum=FALSE )
{
  levs <- if (is.null(levels(cls))) sort(unique(cls)) else levels(cls)
  cls1 <- cls[sub]
  cls1.lev <- levs[match(sort(unique(cls1)),sort(unique(cls)))]
  if (do.renum) cls1 <- match(cls1,unique(cls1))-1
  levels(cls1) <- cls1.lev
  cls1
}
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
change.description <- function( dat, gene2desc, ignore=FALSE )
{
  idx <- match(genenames(dat),gene2desc[,1])
  
  if (any(is.na(idx))) {
    if (ignore) {
      warning(sum(is.na(idx)), " probes w/o description, discarding.\n" )
      idx <- idx[!is.na(idx)]
    }
    else {
      stop( sum(is.na(idx)), " probes w/o description, exiting" )
    }
  }
  newidx <- if (ignore) match(gene2desc[idx,1],rownames(dat@signal)) else 1:nrow(dat)
  newdat <- if(class(dat)=="resdata")
    new("resdata",
        signal=dat@signal[newidx,],
        calls=dat@calls[newidx,],
        description=sapply(gene2desc[idx,2],as.character),
        scale=dat@scale)
  else if (class(dat)=="gctdata")
    new("gctdata",
        signal=dat@signal[newidx,],
        description=sapply(gene2desc[idx,2],as.character))
  else
    stop( "unrecognized data format: ", class(dat) )
  
  newdat
}
xcel2cls <- function( X, levs=NULL, do.sort=FALSE )
{
  if (is.null(levs))
    levs <- unique(X)
  X <- match(X,levs)-1
  ord <- if (do.sort) order(X) else 1:length(X)
  X <- X[ord]
  levels(X) <- levs
  if (do.sort)
    return(list(cls=X,order=ord))
  else
    return(X)
}
factor2cls <- function(FCT)
{
  cls <- match(FCT,levels(FCT))-1
  levels(cls) <- levels(FCT)
  cls
}
## READ GMT
##
## Read gmt file into a named list
##
read.gmt <- function( gmtfile, verbose=TRUE )
{
  gsets <- lapply(scan(gmtfile,what="character",sep="\n",quiet=TRUE),
                  function(z) unlist(strsplit(z,"\t"))[-2])
  names(gsets) <- sapply(gsets,function(z) z[1])
  
  ## *** IMPORTANT: all gene names are 'upper-cased' and replicates are removed ***
  gsets <- lapply(gsets,function(z) {z <- z[-1]; unique(toupper(z[z!=""]))}) # <== upper-case + removal
  gsets
}
## GMT 2 TABLE
##
## create a 0-1 table of genesets from a GSEA '.gmt' file. The output is in
## the format appropriate for hyper.enrichment (the 'categories' argument)
##
gmt2table <- function
(
 gmtfile,
 verbose=TRUE,
 do.save=FALSE,   # save binary format of 0-1 table
 force.read=FALSE # force to read 'gmt' table even if binary already available
)
{
  binfile <- gsub(".gmt",".RData",gmtfile)
  if ( !force.read && file.access(binfile)==0 ) {
    VERBOSE(verbose, "binary object found, loading...\n")
    return( load.var(binfile) )
  }
  # creating a list of genesets
  # *** IMPORTANT: all gene names are 'upper-cased' and replicates are removed ***
  #
  gsets <- read.gmt(gmtfile,verbose=verbose)
  items <- sort(unique(unlist(gsets)))
  pathways <- matrix(0,length(items),length(gsets),
                     dimnames=list(items,names(gsets)))
  percent <- 0.1
  for ( i in 1:length(gsets) )
  {
    gidx <- match(gsets[[i]],rownames(pathways))
    if ( any(is.na(gidx)) )
      stop("something wrong")
    pathways[gidx,i] <- 1
    if (verbose && i>=round(percent*length(gsets)) ) {
      VERBOSE(verbose,round(percent*100),"% ",sep="")
      percent <- percent+0.1
    }
  }
  VERBOSE(verbose,"done, [",paste(dim(pathways),collapse="x"),"] table created.\n",sep="")
  if ( do.save ) {
    VERBOSE(verbose,"Saving binary object .. ")
    save(pathways,file=binfile)
    VERBOSE(verbose,"done.\n")
  }
  return(pathways)
}
