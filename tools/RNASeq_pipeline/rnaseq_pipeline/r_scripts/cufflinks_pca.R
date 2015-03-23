

#Parameter handling
args <- commandArgs(trailingOnly = TRUE)

if (length(args)%%2 !=0){
   print("You need to specify the annotation (-a), the raw counts (-c) and the output directory -o")
   quit(save='no')
}
keys<-args[(1:(length(args)/2))*2-1]
values<-args[(1:(length(args)/2))*2]

#check parameters
chkPars<-function(x,keys,values){
   if (!x %in% keys){
      print(paste('You need to specify the',x,'parameter'))
      quit(save='no')
   }
   return (values[grep(x,keys)])
}

out_dir<-chkPars('-o',keys,values)
counts_file<-chkPars('-c',keys,values)
annot_file<-chkPars('-a',keys,values)

#out_dir='./'
#counts_file<-'cufflinks_counts_fpkm.txt'
#annot_file<-'sample_info.txt'

#read files
annot<-read.table(annot_file,header=T,sep='\t',as.is=T)
rownames(annot)<-gsub('[-\\.]','_',annot$sample_name)

counts<-read.table(counts_file,header=T,sep='\t',as.is=T)
counts<-counts[-1,]
mat<-counts[,-(1:2)]
colnames(mat)<-gsub('[-\\.]','_',colnames(mat))
mat<-apply(mat,2,as.numeric)

#run PCA for all the samples
data<-mat
if ('clickme' %in% rownames(installed.packages())){
   library(clickme)
   data<-data[rowSums(data>100)>=2,]

   # do a PCA
   all<-prcomp(t((data)))
   all<-all$x[,1:2]
   all<-all/apply(all,2,max)

   clickme("points", 
        all[,1],
        all[,2],
        color_groups=annot[,ncol(annot)],
        names = rownames(all),
        title = "Principle Component Analysis",
        xlab = "PCA 1", 
        ylab = "PCA 2",
        file = "pca.html",
        dir=paste(out_dir,'report/cufflinks/',sep=''))$hide()
      
}