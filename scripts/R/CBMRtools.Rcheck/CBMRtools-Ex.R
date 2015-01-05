pkgname <- "CBMRtools"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "CBMRtools-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('CBMRtools')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("heatmap.ggplot2")
### * heatmap.ggplot2

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: heatmap.ggplot2
### Title: heatmap.ggplot2
### Aliases: heatmap.ggplot2

### ** Examples

#Use example data #2, for data set information: ?eSet2
data(eSet2)

#Most basic plot
p1<-heatmap.ggplot2(eSet2, col.clust=FALSE, row.clust=FALSE)

#More advanced plot
p2<-heatmap.ggplot2(eSet=eSet2, col.clust=TRUE, row.clust=TRUE,
col.lab=c("HER2_status", "ER_status", "PR_status", "TN_status"), row.lab="",
heatmap.y.text=FALSE, heatmap.x.text=FALSE,
heatmap.colorlegend.name="RNASeq_expression",
title.text="TCGA BRCA log2 RNA-seq expression, z-score row normalized",
col.legend.name=c("HER2_status", "ER_status", "PR_status", "TN_status"),
row.legend.name="",
row.scaling="z-score.capped",
z.norm=FALSE,
cuttree.col=4, cuttree.row=6,
verbose=FALSE, show=FALSE)

#Display and saving options
print(p2) #to display in viewer if working in RStudio
ggsave(p2, file="example_plot.pdf")
ggsave(p2, file="example_plot.jpg")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("heatmap.ggplot2", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("variation.filter")
### * variation.filter

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: variation.filter
### Title: variation.filter
### Aliases: variation.filter

### ** Examples

# select the top 100 genes by median absolute deviation (mad) and also generate plot

#Use example data #2, for data set information: ?eSet2
data(eSet2)

png('mydata.filtered.png')
variation.filter(eset2,score='mad',ngenes=100,do.plot=TRUE)
dev.off()
...



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("variation.filter", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
