#' @export

########### MAIN #########
##########################

if ( !exists("CBMDEV") ) {
  CBMDEV <- Sys.getenv('CBMDEV')
  if (CBMDEV=="") stop( "Use 'setenv CBMDEV ..' to set CBMrepository's base directory" )
}
if ( !exists("CBMMLAB") ) {
  CBMMLAB <- Sys.getenv('CBMMLAB')
  if (CBMMLAB=="") stop( "Use 'setenv CBMMLAB ..' to set CBMrepository's base directory" )
}
source(paste(CBMDEV,"R/misc.R",sep="/"))
source(paste(CBMDEV,"R/misc.math.R",sep="/"))
source(paste(CBMDEV,"R/broad.file.formats.R",sep="/"))
source(paste(CBMDEV,"R/permute.array.R",sep="/"))
source(paste(CBMDEV,"R/perm.1side.R",sep="/"))
source(paste(CBMDEV,"R/perm.2side.R",sep="/"))
source(paste(CBMDEV,"R/diffanal.scores.R",sep="/"))
source(paste(CBMDEV,"R/nn.analysis.R",sep="/"))
