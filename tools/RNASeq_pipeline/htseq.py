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

import module_helper, helper, subprocess, os

def init(param):
   module_helper.checkParameter(param,key='HTSeq_exec',dType=str)
   module_helper.checkParameter(param,key='sam_exec',dType=str)
   module_helper.checkParameter(param,key='HTSeq_s',dType=str)
   module_helper.checkParameter(param,key='HTSeq_t',dType=str)
   module_helper.checkParameter(param,key='HTSeq_m',dType=str,allowed=['union', 'intersection-strict','intersection-nonempty'])
   module_helper.checkParameter(param,key='HTSeq_id',dType=str)
   module_helper.checkParameter(param,key='HTSeq_gft',dType=str,checkFile=True)
   module_helper.checkParameter(param,key='Rscript_exec',dType=str)
   
def createESet(count_file,pheno_file,param):
    #create a Bioconductor ExpresionSet
    cmd=[param['Rscript_exec'],'createRawCountESet.R','-c',count_file,'-a',pheno_file,'-o',param['working_dir']]
    output,error = subprocess.Popen(cmd ,stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()
    helper.writeLog('Creating ESet ... \n',param)
    

def copyAndReadStatFiles(param):
    #if there is no htseq directory in the report make one
    htseq_dir=param['working_dir']+'report/htseq/'
    if not os.path.exists(htseq_dir):
        os.makedirs(htseq_dir)
 
    #get the files that are actually in the output directory       
    call=['cp','-R',param['working_dir']+'results/htseq/',param['working_dir']+'report/']
    output,error = subprocess.Popen(call,stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()
    print error
    print output
    
    htseq_file=param['working_dir']+'results/htseq/htseq_stats.txt'    
    #extract table
    table=[]
    fh=open(htseq_file)
    #header
    table.append(fh.readlines()[0].rstrip().split('\t'))    
    table[0]=table[0][1:]
    fh.close()
   
    #total number of aligned reads
    tot_reads=['Total number of aligned reads']+param['bam_qc']['total_aligned_reads']
    table.append(tot_reads)    

    fh=open(htseq_file)
    for line in (fh.readlines()[1:]):
       cur_line=line.rstrip().split('\t')
       perc=[cur_line[0]]+module_helper.getPercentage(cur_line[1:],tot_reads[1:],len(cur_line)-1)   
       table.append(perc)
              
    fh.close()
    return table
 


def report(param):
    #report only if there were actually results
    if os.path.exists(param['working_dir']+'deliverables/htseq_raw_counts.txt'):
        param['report'].write('<center><br><br><br><br><h2>HTSeq statistics</h2>')    
        table=copyAndReadStatFiles(param)
        module_helper.writeHTMLtable(param,table,out=param['report'],cell_width=80,fCol_width=150,deg=315)
        param['report'].write('<a href="htseq/htseq_stats.txt">HTSeq statistics as tab delimited txt file</a>')
        


def finalize(param,input_files='count_files'):
    helper.writeLog('Collecting HTSeq raw counts ... \n',param)
    #extracts the counts from the htseq output
    import csv
 
    #check which of these files are actually available
    working_files = [iFile for iFile in param[input_files] if iFile !='']
    
    if len(working_files)>0:
        #get gene annotation
        csv_file = open(working_files[0])
        csv_reader = csv.reader(csv_file, delimiter='\t')
        counts = [row[0] for row in csv_reader]
        csv_file.close()
    
        #get all the expression values
        header='ENS_ID'
        for idx in range(param['num_samples']):
            if param[input_files]!='':            
                header=header+'\t'+param['stub'][idx]
                csv_file = open(param[input_files][idx])
                csv_reader = csv.reader(csv_file, delimiter='\t')
                i=0
                for row in csv_reader:
                    counts[i] = counts[i]+'\t'+row[1]
                    i+=1 
                csv_file.close()
    
        #output the file
        out_file=param['working_dir']+'deliverables/htseq_raw_counts.txt'
        out_handle=open(out_file,'w')
        out_handle.write(header+'\n')
        #drop the last 5 lines since these provide only a summary statistic
        for i in range(len(counts)-5): out_handle.write(counts[i]+'\n')
        out_handle.close()
        
        #output_phenotype_file
        helper.writeLog('Writing phenotype data ... \n',param)
        module_helper.output_sample_info(param)
        
        #create an eSet:
        createESet(out_file,param['pheno_file'],param)
        
        #write summary stats        
        out_handle=open(param['working_dir']+'results/htseq/htseq_stats.txt','w')
        out_handle.write(header+'\n')
        for i in range(len(counts)-5,len(counts)): out_handle.write(counts[i]+'\n')
        out_handle.close()
        
    else:
        print('HTseq was not run successfully on any of the files..\n')
        
        
        
if __name__ == "__main__":
    import subprocess, sys, os
    param=module_helper.initialize_module()

    #build htseq-count call:
    call1 = [param['sam_exec'],'view',param['working_file']]
    call2 = [param['HTSeq_exec'],'-s',param['HTSeq_s'],'-t',param['HTSeq_t'],'-i',param['HTSeq_id'],'-m',param['HTSeq_m'],'-',param['HTSeq_gft']]
  
    #function calls
    param['file_handle'].write('Pipe CALL 1: '+' '.join(call1)+'\n')
    param['file_handle'].write('Pipe CALL 2: '+' '.join(call2)+'\n')
    
    p1 = subprocess.Popen(call1, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(call2, stdin=p1.stdout, stdout=subprocess.PIPE)
    p1.stdout.close()
    output,error = p2.communicate()
        
    #error handling
    if output=='':
        param['file_handle'].write(error+'\n') 
        sys.exit(0)
    
    #write output
    outfile=param['module_dir']+param['stub'][param['file_index']]+'.txt'
    handle=open(outfile,'w')
    handle.write(output)
    handle.close()
    
    #wrap up and return the current workingfile
    module_helper.wrapup_module(param,[outfile]) 
     


