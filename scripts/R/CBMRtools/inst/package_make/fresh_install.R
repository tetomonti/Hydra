#installing dependencies not from cran
install.packages("devtools")
require(devtools)

install_github('hadley/staticdocs')
install_github("andrie/ggdendro")

#installing CBMRtools
package.dir <- normalizePath("../../../CBMRtools")
install.packages(package.dir, type = "source", repo = NULL) 

