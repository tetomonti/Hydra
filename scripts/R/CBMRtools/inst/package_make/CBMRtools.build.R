## Sequence of commands to check and build the package 
## and generate html page documents

#install.packages(devtools)
#install_github('hadley/staticdocs')

require(devtools)
require(staticdocs)

package.dir <- normalizePath("../../../CBMRtools")
document(package.dir) #creates help pages
check(package.dir) #package checking
load_all(package.dir) #building package
library(CBMRtools)

#generate html pages
setwd(package.dir)
build_site(pkg = package.dir, examples = TRUE, launch = TRUE)