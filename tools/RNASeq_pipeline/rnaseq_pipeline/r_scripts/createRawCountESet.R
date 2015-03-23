

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

#read files
annot<-read.table(annot_file,header=T,sep='\t',as.is=T)
rownames(annot)<-gsub('[-\\.]','_',annot$sample_name)

counts<-read.table(counts_file,header=T,sep='\t',as.is=T)
rownames(counts)<-counts[,1]
counts<-counts[,-1]
colnames(counts)<-gsub('[-\\.]','_',colnames(counts))


#create and save eSet
require(Biobase)
metadata<-data.frame(labelDescription=colnames(annot),row.names=colnames(annot))               
phenoData<-new("AnnotatedDataFrame", data=annot, varMetadata=metadata)   
expr.data<-new("ExpressionSet", 
               exprs=as.matrix(counts), 
               phenoData=phenoData, 
               annotation='RNASeq raw counts')
saveRDS(expr.data,file=paste(out_dir,'deliverables/rawCounts.RDS',sep=''))

