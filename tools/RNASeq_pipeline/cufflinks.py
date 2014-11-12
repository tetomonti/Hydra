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

import module_helper, helper, os, subprocess

def init(param):
   module_helper.checkParameter(param,key='cufflinks_exec',dType=str)
   module_helper.checkParameter(param,key='cufflinks_G',dType=str,checkFile=True)
   module_helper.checkParameter(param,key='cufflinks_b',dType=str,checkFile=True)
   module_helper.checkParameter(param,key='cufflinks_library_type',dType=str)
   module_helper.checkParameter(param,key='cufflinks_compatible_hits',dType=str)
   module_helper.checkParameter(param,key='cufflinks_total_hits',dType=str)
   module_helper.checkParameter(param,key='cufflinks_N',dType=str)
   module_helper.checkParameter(param,key='cufflinks_u',dType=str)   
   
   

def report(param):
    param['report'].write('<center><br><br><h2>Cufflinks results</h2>')    

    #make cufflinks directory
    cufflinks_dir=param['working_dir']+'report/cufflinks/'
    if not os.path.exists(cufflinks_dir):
        os.makedirs(cufflinks_dir)

    #output_phenotype_file
    module_helper.output_sample_info(param)
        
    #run R script that creates a PCA
    counts=param['working_dir']+'deliverables/cufflinks_counts_fpkm.txt'
    cmd=[param['Rscript_exec'],'cufflinks_pca.R','-c',counts,'-a',param['pheno_file'],'-o',param['working_dir']]
    output,error = subprocess.Popen(cmd ,stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()
    param['report'].write('<a href="cufflinks/pca.html">PCA</a>')
        
    
def finalize(param,input_files='count_files'):
    helper.writeLog('Collecting Cufflinks FPKM data ... \n',param)
    #extracts the counts from the cufflinks output
    import csv
 
    #check which of these files are actually available
    working_files = [iFile for iFile in param[input_files] if iFile !='']
    
    if len(working_files)>0:
        #get gene annotation
        csv_file = open(working_files[0])
        csv_reader = csv.reader(csv_file, delimiter='\t')
        counts = [row[0]+'\t'+row[4] for row in csv_reader]
        csv_file.close()
    
        #get all the expression values
        header='ENS_ID\tSymbol'
        for idx in range(param['num_samples']):
            if param[input_files]!='':            
                header=header+'\t'+param['stub'][idx]
                csv_file = open(param[input_files][idx])
                csv_reader = csv.reader(csv_file, delimiter='\t')
                i=0
                for row in csv_reader:
                    counts[i] = counts[i]+'\t'+row[9]
                    i+=1 
                csv_file.close()
    
        #output the file
        out_handle=open(param['working_dir']+'deliverables/cufflinks_counts_fpkm.txt','w')
        out_handle.write(header+'\n')
        for i in range(len(counts)): out_handle.write(counts[i]+'\n')
        out_handle.close()
    else:
        print('Cufflinks was not run successfully on any of the files..\n')
    
    
if __name__ == "__main__":
    import subprocess, sys, os
    param=module_helper.initialize_module()

    #run create output directory
    outdir = param['module_dir']+param['stub'][param['file_index']]+'/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    #build cufflinks call:
    call = param['cufflinks_exec']+' '
    
    #optional parameters
    if param['cufflinks_compatible_hits']=='active':
        call=call+'--compatible-hits-norm '
    if param['cufflinks_N']=='active':
        call=call+'-N '
    if param['cufflinks_u']=='active':
        call=call+'-u '
    if param['cufflinks_total_hits']=='active':
        call=call+'--total-hits-norm '
        
    #finish call
    call=call+' --no-update-check --library-type '+param['cufflinks_library_type']+' -o '+outdir+' -p '+param['num_processors']+' -G '+param['cufflinks_G']+' '+ param['working_file']
    
    param['file_handle'].write('CALL: '+call+'\n')
    output,error = subprocess.Popen(call.split() ,stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()
    param['file_handle'].write(output)
    
    #printing the error blows the log file up to several megabytes
    #param['file_handle'].write(error)
     
    #error handling
    if not os.path.exists(outdir+'genes.fpkm_tracking'):
        param['file_handle'].write('Error there was not cuffinks output\n') 
        sys.exit(0)
    
    #wrap up and return the current workingfile
    module_helper.wrapup_module(param,[outdir+'genes.fpkm_tracking']) 


