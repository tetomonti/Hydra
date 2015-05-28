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

#extract the relevant parameters
out_dir<-chkPars('-o',keys,values)
counts_file<-chkPars('-c',keys,values)
annot_file<-chkPars('-a',keys,values)
paired<-chkPars('-p',keys,values)
stub<-chkPars('-s',keys,values)


#actual values for testing purposes
#out_dir='./'
#counts_file<-'deliverables/cufflinks_counts_fpkm.txt'
#annot_file<-'deliverables/sample_info.txt'
#stub='cufflinks'
#paired='FALSE'


#read phenotype file
annot<-read.table(annot_file,header=T,sep='\t',as.is=T)
rownames(annot)<-gsub('[-\\.]','_',annot$sample_name)


#read the raw counts file
counts<-read.table(counts_file,header=T,sep='\t',as.is=T)
counts<-counts[-1,]
mat<-as.matrix(counts[,-(1:2)])
colnames(mat)<-gsub('[-\\.]','_',colnames(mat))
mat<-apply(mat,2,as.numeric)
data<-mat

###################################################################
#PCA and clickme
###################################################################

#function to create a 
plotCov<-function(idx,annot,all,stub,con){
   current_filename <- paste0(stub,'_',colnames(annot)[idx])
   code<-clickme("points", 
                 all[,1],
                 all[,2],
                 color_groups=annot[,idx],
                 names = rownames(all),
                 title = colnames(annot)[idx],
                 xlab = "PCA 1", 
                 ylab = "PCA 2",
                 file_path = paste0(out_dir,'report/clickme/',current_filename,'.html'))
   write(paste0('<iframe width="1100" height="850" src="',
               '../clickme/',
               current_filename,
               '.html" frameborder=0> </iframe>'),
         con)
}

#run PCA for all the samples
#only if clickme is installed and there is more than 1 sample
if ('clickme' %in% rownames(installed.packages()) & ncol(data)>1){
   library(clickme)
   #if the clickme directory does not work create it
   dir.create(file.path(paste0(out_dir,'report/'), 'clickme'), showWarnings = FALSE)
      
   #remove the raw file names
   if (paired == 'TRUE'){
      annot<-as.matrix(annot[,-(1:2)])
   }else{
      annot<-as.matrix(annot[,-1])
   }
   
   # do a PCA
   all<-prcomp(t(log2(data+1)))
   all<-all$x[,1:2]
   all<-all/apply(all,2,max)

   #plot the pca   
   con<-file(paste0(out_dir,'report/',stub,'/',stub,'_pca.html'),open='w')
   
   #always plot the sample names, just in case
   plotCov(1,annot,all,stub,con)
   annot<-annot[,-1]
   
   #and then plot all covariates that have more than 1 and equal or less than 10 levels/classes
   no_levels<-apply(annot,2,function(x)length(unique(x)))
   valid_indices<-(1:ncol(annot))[no_levels>1 & no_levels<=10]
   
   #only if there is actually something to plot
   if (length(valid_indices)>0){   
      sapply(valid_indices,plotCov,annot,all,stub,con)
   }
   
   #close the html file
   close(con) 
}


###################################################################
#Boxplots
###################################################################

exp<-as.matrix(log2(data+1))
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
