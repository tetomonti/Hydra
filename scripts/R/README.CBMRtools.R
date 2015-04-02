CBMGIT <- Sys.getenv('CBMGIT')
if (CBMGIT=="") stop( "Use 'setenv CBMGIT ..' to set CBMgithub's base directory" )

setwd(paste(CBMGIT,'/scripts/R/CBMRtools/inst/package_make',sep=''))
source(“CBMRtools.build.R”) 

