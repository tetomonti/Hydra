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

import module_helper
import helper
import fastqc
import trimmer
import tophat
import bamqc
import cufflinks 
import matched_pairs
import cutadapt
import htseq
import sys

    
def initialize_all(param):
    # this function calls the initialize functions of every module
    trimmer.init(param)
    cutadapt.init(param)
    matched_pairs.init(param)
    fastqc.init(param)
    tophat.init(param)
    bamqc.init(param)
    cufflinks.init(param)  
    htseq.init(param)  

    
def run_all(param):
    
    #is the input bam or 
    if param['aligner']=='skip':
       param['bam_files']=param['raw_files'][:]
    else:
        #preprocessing fastq file
        helper.submit_job(param, 'fastqc.py',input_files='raw_files')  
        if param['do_trimming']:
            helper.submit_job(param,'trimmer.py',input_files='raw_files',output_files='fastq_files') 
        else:
           param['fastq_files']=param['raw_files'][:]
           if param['paired']:
               param['fastq_files2']=param['raw_files2'][:]
        helper.submit_job(param, 'cutadapt.py',    input_files='raw_files',output_files='fastq_files',environment='modules_python2.7.5')
        if param['paired']:
               helper.submit_job(param, 'matched_pairs.py', input_files='fastq_files',output_files='fastq_files',environment='modules_python2.7') 
        helper.submit_job(param, 'fastqc.py',    input_files='fastq_files') 
        
        #do alignment if it's not just a fastqc run
        if not param['QC_and_trim_only']:
            if param['aligner']=='tophat':
                #running the aligner
                helper.submit_job(param, 'tophat.py',    input_files='fastq_files', output_files='bam_files', cores=param['qsub_num_processors'])  
            else:
                helper.writeLog('The selected aligner does not exist.',param)
                sys.exit(0)
            
    if not param['QC_and_trim_only']:
        #Bamqc
        helper.submit_job(param, 'bamqc.py',     input_files='bam_files',   environment='modules_python2.7.5') 
        
        #Getting the counts:
		#Added code (conditional statements)
        if param['run_cufflinks']:
          helper.submit_job(param, 'cufflinks.py', input_files='bam_files',   output_files='count_files')  
          cufflinks.finalize(param,input_files='count_files')
        if param['run_htseq']:
          helper.submit_job(param, 'htseq.py', input_files='bam_files',   output_files='count_files', environment='modules_python2.7')  
          htseq.finalize(param,input_files='count_files')


def report_all(param):
    if param['aligner']!='skip':
        fastqc.report(param,input_files='raw_files',header='FastQC results on the raw data')
        fastqc.report(param,input_files='fastq_files',header='FastQC results after preprocessing')
    bamqc.report(param)
	#Added code (conditional statements)
    if param['run_cufflinks']:
        cufflinks.report(param)
    if param['run_htseq']:
        htseq.report(param)

    
    
    
