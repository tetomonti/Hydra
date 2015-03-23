import rnaseq_pipeline.module_helper
module_helper = rnaseq_pipeline.module_helper
import os

def init(param):
   module_helper.checkParameter(param,key='trimming_exec',dType=str)
   module_helper.checkParameter(param,key='trimming_mode',dType=str,allowed=['-l','-f','-t'])
   module_helper.checkParameter(param,key='trimming_value',dType=str)
   module_helper.checkParameter(param,key='trimming_quality_value',dType=str,optional=True)
   module_helper.checkParameter(param,key='do_trimming',dType=bool)
   
   
   

def run_trimmer(param,infile,outfile):

    #unzip the fastq file
    temp_file = "'"+param['module_dir']+param['stub'][param['file_index']]+".fastq'"
    call = 'gunzip -c ' + param[infile] + ' >'+temp_file
    retvalue = os.system(call)

    #run the trimmer
    param['file_handle'].write('CALL: '+call+'\n')
    call = param['trimming_exec']+' '+ param['trimming_mode'] + ' ' + param['trimming_value'] + ' ' + param['trimming_quality_value'] +' -z -i '+ temp_file +' -o '+param[outfile]
    output,error = subprocess.Popen(call ,shell=True, stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()
    
    #remove temporary file
    os.remove(temp_file)
        
    param['file_handle'].write(error)
    param['file_handle'].write(output)
 
    # Error handling
    if (not os.path.exists(param[outfile])) or os.stat(param[outfile]).st_size<1000000:
        sys.exit(0)


def main():
    import subprocess, sys, os
    param=module_helper.initialize_module()

    #run fastqc
    param['newfile']=param['module_dir']+param['stub'][param['file_index']]+'.trimmed.fastq.gz'
    run_trimmer(param,'working_file','newfile')
    
 
    if not param['paired']:
        module_helper.wrapup_module(param,[param['newfile']])
    #calling it on the second fastq file if it is paired    
    else:    
        param['newfile2']=param['module_dir']+param['stub'][param['file_index']]+'.trimmed.2.fastq.gz'
        run_trimmer(param,'working_file2','newfile2')
        module_helper.wrapup_module(param,[param['newfile'],param['newfile2']]) 
        
