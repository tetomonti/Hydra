import rnaseq_pipeline.module_helper
module_helper = rnaseq_pipeline.module_helper
import os
import sys
import subprocess

def init(param):
   module_helper.checkParameter(param,key='cutadapt_exec',dType=str)
   module_helper.checkParameter(param,key='cutadapt_first_adapter',dType=str)
   if param['paired']:
       module_helper.checkParameter(param,key='cutadapt_second_adapter',dType=str)
   module_helper.checkParameter(param,key='cutadapt_m',dType=str)
   
   
def run_cutadapt(param,infile,outfile,adapter):
    #run the clipper

    call = param['cutadapt_exec']+' -m '+param['cutadapt_m']+' -a '+adapter+' -o "'+outfile+'" "'+param[infile]+'"'
    print '\n'
    print call
    print '\n'

    output,error = subprocess.Popen(call ,stdout = subprocess.PIPE, shell=True, stderr= subprocess.PIPE).communicate()
    print output
    print error   

    param['file_handle'].write(error)
    param['file_handle'].write(output)

    with open(outfile+'.txt','w') as f:
        f.write(output)

    # Error handling
    if (not os.path.exists(outfile)) or os.stat(outfile).st_size<10000:
        sys.exit(0)


def main():
    import subprocess, sys, os
    param=module_helper.initialize_module()

    #run fastqc
    outfile = param['module_dir']+param['stub'][param['file_index']]+'.clipped.fastq.gz'
    run_cutadapt(param,'working_file',outfile,param['cutadapt_first_adapter'])
     
    if not param['paired']:
        module_helper.wrapup_module(param,[outfile])
    #calling it on the second fastq file if it is paired    
    else:    
        outfile2 = param['module_dir']+param['stub'][param['file_index']]+'.clipped.2.fastq.gz'
        run_cutadapt(param,'working_file2',outfile2,param['cutadapt_second_adapter'])
        module_helper.wrapup_module(param,[outfile,outfile2]) 
