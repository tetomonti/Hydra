#Copyright 2015 Daniel Gusenleitner, Stefano Monti

#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at

#    http://www.apache.org/licenses/LICENSE-2.0

#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.


############
#these are test parameters
#out_dir='./'
#counts_file<-'deliverables/htseq_raw_counts.txt'
#annot_file<-'deliverables/sample_info.txt'
#stub='htseq'
#paired=TRUE


###################################################################
#Parameter handling
###################################################################
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 10){
   print("You need to specify the annotation (-a), the raw counts (-c),")
   print("the working directory (-o), whether the reads are paired (-p), and the stub name (-s)")
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
paired<-chkPars('-p',keys,values)
stub<-chkPars('-s',keys,values)


###################################################################
#Output raw files
###################################################################

#read files
annot<-read.table(annot_file,header=T,sep='\t',as.is=T)
rownames(annot)<-gsub('[-\\.]','_',annot$sample_name)

counts<-read.table(counts_file,header=T,sep='\t',as.is=T)
gene_names<-counts[,1]
sample_names<-colnames(counts)[-1]
counts<-as.matrix(counts[,-1])
colnames(counts)<-gsub('[-\\.]','_',sample_names)
rownames(counts)<-gene_names

#get annotation only for samples that were sucessfully processed
if (paired){
   idx<-3
}else{
   idx<-2
}
#fix potential sample name issue that arises when a sample starts with a number or an underscore
#R prepends a 'X' in this case
names1<-gsub('[/_-]','.',annot[,idx])
odd_ids<-c(grep('^\\.',names1),grep('^\\d',names1))
if (length(odd_ids)>0){
   names1[odd_ids]<-paste0('X',names1[odd_ids])
}
names2<-gsub('[/_-]','.',colnames(counts))
annot<-annot[match(names1,names2),]
colnames(counts)<-rownames(annot)


#create and save eSet
suppressMessages(require(Biobase))
metadata<-data.frame(labelDescription=colnames(annot),row.names=colnames(annot))               
phenoData<-new("AnnotatedDataFrame", data=annot, varMetadata=metadata)   
expr.data<-new("ExpressionSet", 
               exprs=as.matrix(counts), 
               phenoData=phenoData, 
               annotation='RNASeq raw counts')
saveRDS(expr.data,file=paste0(out_dir,'deliverables/',stub,'_raw_counts.RDS'))

###################################################################
#Normalizing and output counts
###################################################################
suppressMessages(require(edgeR))

eSet<-expr.data
exprs(eSet)<-cpm(eSet)
saveRDS(eSet,file=paste0(out_dir,'deliverables/',stub,'_normalized_counts.RDS'))


###################################################################
#PCA and clickme
###################################################################

#function to create a 
plotCov<-function(idx,variances,annot,all,stub,con){
   current_filename <- paste0(stub,'_',colnames(annot)[idx])
   code<-clickme("points", 
                 all[,1],
                 all[,2],
                 color_groups=annot[,idx],
                 names = rownames(all),
                 title = colnames(annot)[idx],
                 x_title = paste0("PCA 1 (",variances[1],"%)"), 
                 y_title = paste0("PCA 2 (",variances[2],"%)"),
                 file_path = paste0(out_dir,'report/clickme/',current_filename,'.html'))
   write(paste0('<iframe width="1100" height="850" src="',
                '../clickme/',
                current_filename,
                '.html" frameborder=0> </iframe>'),
         con)
}


#run PCA for all the samples
#only if clickme is installed and there is more than 1 sample
if ('clickme' %in% rownames(installed.packages()) & ncol(eSet)>1){
  suppressMessages(library(clickme))
  #if the clickme directory does not work create it
  dir.create(file.path(paste0(out_dir,'report/'), 'clickme'), showWarnings = FALSE)
  
  header <- colnames(annot)
  #remove the raw file names
  if (paired == 'TRUE'){
    annot<-as.matrix(annot[,-(1:2)])
    colnames(annot)<-header[-(1:2)]
  }else{
    annot<-as.matrix(annot[,-1])
    colnames(annot)<-header[-1]
  }
  
  # do a PCA
  all<-prcomp(t(exprs(eSet)))
  y<-cov(all$x)
  variances<-c(y[1,1]/sum(diag(y)),y[2,2]/sum(diag(y))) * 100
  variances<-round(variances,digits=2)
  all<-all$x[,1:2]
  all<-all/apply(all,2,max)
  
  #plot the pca   
  con<-file(paste0(out_dir,'report/',stub,'/',stub,'_pca.html'),open='w')
  
  #and then plot all covariates that have more than 1 and equal or less than 10 levels/classes
  no_levels<-apply(annot,2,function(x)length(unique(x)))

  valid_indices<-(1:ncol(annot))[no_levels>1 & no_levels<=10]
  if (length(valid_indices) == 0){
     valid_indices <- 1
  }  

  #only if there is actually something to plot
  sapply(valid_indices,plotCov,variances,annot,all,stub,con)
  
  #close the html file
  close(con) 
}



###################################################################
#Additional QC
###################################################################

exp<-as.matrix(exprs(eSet))
exp<-as.matrix(exp[order(apply(exp,1,mad),decreasing = T)[1:5000],])

png(paste0(out_dir,'report/',stub,'/',stub,'_boxplot.png'),
    height=800,
    width=1024+10*ncol(exp))
par(mar=c(10, 4.1, 4.1, 2.1))
labels <-colnames(exp)
boxplot(exp,col="lightgray",xaxt="n",xlab = "")
axis(1, labels = FALSE)
text(x =  seq_along(labels), y = par("usr")[3] - 1, srt = 90, adj = 1,
     labels = labels, xpd = TRUE)
invisible(dev.off())







