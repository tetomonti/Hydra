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
import featureCount
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
    featureCount.init(param)  

    
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
#Added code
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
        #Added code (conditional statements to accomodate multiple read quantifiers)
        if param['run_cufflinks']:
          helper.submit_job(param, 'cufflinks.py', input_files='bam_files',   output_files='count_files')  
          cufflinks.finalize(param,input_files='count_files')
        
        if param['run_htseq']:
	  helper.submit_job(param, 'htseq.py', input_files='bam_files',   output_files='count_files', environment='modules_python2.7')  
          htseq.finalize(param,input_files='count_files')

        if param['run_featureCount']:		  
          helper.submit_job(param,'featureCount.py',input_files='bam_files', output_files='count_files',environment='modules_python2.7')
	  featureCount.finalize(param,input_files='count_files')

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
    if param['run_featureCount']:
     featureCount.report(param)
    
    
    
