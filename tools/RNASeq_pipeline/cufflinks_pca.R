#  Copyright (c) 2014, Boston University. All rights reserved.
#  
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are met: 
#  
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer. 
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution. 
#  
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
#  ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
#  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#  
#  The views and conclusions contained in the software and documentation are those
#  of the authors and should not be interpreted as representing official policies, 
#  either expressed or implied, of Boston University.
#  
#  Authors:
#    Daniel Gusenleitner [1,2], Vinay Kartha [1,2], Francesca Mulas [2], 
#    Yuxiang Tan [1,2], Liye Zhang [2], Stefano Monti [1,2]
#
#  [1] Bioinformatics Program, Boston University
#  [2] Center for Computational Biomedicine, Boston University  
#  



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