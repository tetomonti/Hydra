## Sequence of commands to check and build the package 
## and generate html page documents

#install.packages(devtools)
#require(devtools)
#install_github('hadley/staticdocs')

require(devtools)
require(staticdocs)

package.dir <- normalizePath("../../../CBMRtools")
document(package.dir) # creates help pages
check(package.dir)    # package checking
load_all(package.dir) # building package
library(CBMRtools)

## the directory 'staticdocs' must exist under CBMRtools/inst/ for the
## ..help pages to be installed at CBMRtools/inst/web

STATICDIR <- '../staticdocs'
if ( is.na(file.info(STATICDIR)$isdir) ) system(paste('mkdir',STATICDIR))

## generate html pages
setwd(package.dir)
build_site(pkg = package.dir, examples = TRUE, launch = TRUE)

## install CBMRtools locally
cat("installing package...\n")
install.packages(package.dir, repos = NULL, type = "source")
