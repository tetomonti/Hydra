## Sequence of commands to check and build the package

require(devtools)
package.dir <- "~/Research/CBMgithub/scripts/R/CBMRtools"
document(package.dir) #creates help pages
check(package.dir) #package checking
load_all(package.dir) #building package
library(CBMRtools)


