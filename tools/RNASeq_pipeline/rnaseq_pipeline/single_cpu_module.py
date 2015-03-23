import os, sys,  time, shlex, subprocess
from rnaseq_pipeline import module_helper

def initialize_single_cpu(param):
    #split module list and add module load to the list
    module_helper.checkParameter(param,key='modules',dType=str)
    param['modules']=[]

        
def run_single_job(index, param, job_id, py_file, cores):
    #fetch modules that need to be loaded and add the python command for a single sample
    command_list = param['modules'][:]
    command_list.append(py_file +' -i ' + str(index) + ' -n $NSLOTS')

    #call qsub script
    for cmd in command_list:
        args = shlex.split(cmd)
        output,error = subprocess.Popen(args,stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()

