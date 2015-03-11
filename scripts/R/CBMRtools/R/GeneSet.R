#  Copyright (c) 2012, 2013, Boston University. All rights reserved.
#  
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are met: 
#  
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer. 
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution. 
#  
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
#  ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
#  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#  
#  The views and conclusions contained in the software and documentation are those
#  of the authors and should not be interpreted as representing official policies, 
#  either expressed or implied, of Boston University.
#  
#  Authors:
#    Arjan van der Velde [1], Adam Labadorf [1], Heather Selby [1],
#    Daniel Gusenleitner [1,2], Stefano Monti [1,2]
#  
#  [1] Bioinformatics Program, Boston University
#  [2] Center for Computational Biomedicine, Boston University  
#  

# $LastChangedDate: 2013-10-14 11:36:59 -0400 (Mon, 14 Oct 2013) $
# $LastChangedRevision: 355 $


# GeneSet class. Represents a gene set file/rdata/package data object.
GeneSet <- setClass("GeneSet", 
                    representation(source.file="character",
                                   geneset="list", 
                                   type="character",
                                   verbose="logical", 
                                   do.save="logical", 
                                   name="character"), 
                    prototype(source.file=character(0),
                              geneset=list(), 
                              type="hgnc_symbol",  
                              do.save=T, 
                              verbose=F, 
                              name=NULL), 
                    sealed=FALSE)

# GeneSet constructor.
setMethod("initialize", "GeneSet",
  function(.Object, source.file=NULL, geneset=NULL, type="hgnc_symbol", do.save=T, verbose=F, name=NULL) {
    if(is(geneset, "GeneSet")){
       .Object <- GeneSet(geneset@source.file,
                          geneset@geneset, 
                          type=geneset@type,
                          do.save=geneset@do.save|do.save, 
                          verbose=geneset@verbose|verbose,
                          name=geneset@name)
       
    }else if (is(source.file, "character") && length(source.file) == 1) {
      .Object@verbose <- verbose
      .Object@do.save <- do.save
       
      #checking file format
      file.type<-strsplit(source.file,'\\.')[[1]]
      file.type<-file.type[length(file.type)]
      
      if(file.type=="gmt") {
         if (file.access(source.file)!=0) {
            stop(sprintf("Cannot read gene set file %s", source.file))
         } else {
            x<-read.gmt(source.file,verbose=verbose)
            .Object@geneset <- x
            .Object@name<-source.file
         }         
      }else{
          stop("Only .gmt files are supported.")
      }
   }else{
       stop("Character string or GeneSet expected.")
    }
    return(.Object)
  }
, sealed=FALSE)

# GeneSet: show
setAs("GeneSet", "character", function(from) sprintf("%s (%s): %s", from@name, from@type, from@geneset))
setMethod("show", "GeneSet", function(object) 
   cat('GeneSet Object containing: ',length(object@geneset),'gene sets\n'), 
          sealed=FALSE)
setMethod("print", "GeneSet", function(x, ...) show(x), sealed=FALSE)

# GeneSet: geneSetName
setGeneric("geneSetName", function(object) standardGeneric("geneSetName"))
setGeneric("geneSetName<-", function(object,value) standardGeneric("geneSetName<-"))
setMethod("geneSetName", "GeneSet", function(object) object@name, sealed=FALSE)
setReplaceMethod("geneSetName",
                 "GeneSet",
                 function(object, value) {
                   object@name <- value
                   return(object)
                 },
                 sealed=FALSE)

# GeneSet: Get or set geneset
setGeneric("geneSet", function(object) standardGeneric("geneSet"))
setGeneric("geneSet<-", function(object,value) standardGeneric("geneSet<-"))
setMethod("geneSet", "GeneSet", function(object) return(object@geneset), sealed=FALSE)
setReplaceMethod("geneSet",
                 "list",
                 function(object, value) {
                   object@geneset <- value
                   return(object)
                 },
                 sealed=FALSE)

