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

"""bowtie2 module
This module contains functions for initializing all bowtie2 specific variables and
a wrapper that that runs bowtie2 using those parameters on a single sample.
"""

from hydra_pkg import module_helper as MODULE_HELPER

def init(param):
    """Initialization function that checks the all relevant bowtie2 parameters
    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    
    MODULE_HELPER.check_parameter(param, key='bowtie2_exec', dtype=str)
    MODULE_HELPER.check_parameter(param, key='bowtie2_index', dtype=str)
    MODULE_HELPER.check_parameter(param, key='bowtie2_type', dtype=str)
    MODULE_HELPER.check_parameter(param, key='bowtie2_N', dtype=str)
    MODULE_HELPER.check_parameter(param, key='bowtie2_D', dtype=str)
    MODULE_HELPER.check_parameter(param, key='bowtie2_R', dtype=str)
    MODULE_HELPER.check_parameter(param, key='bowtie2_L', dtype=str)
    MODULE_HELPER.check_parameter(param, key='bowtie2_i', dtype=str)
    MODULE_HELPER.check_parameter(param, key='bowtie2_rdg', dtype=str)
    MODULE_HELPER.check_parameter(param, key='bowtie2_rfg', dtype=str)


def main():
    """Main function that is run on each samples, which in turn calls runs
    bowtie2 on a sample.
    """
    import subprocess
    import sys
    import os
    param = MODULE_HELPER.initialize_module()
    #run create output directory
    outdir = param['module_dir']+param['outstub']+'/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    #build bowtie call:
    call = [param['bowtie2_exec']]
    
    #options
    call.append(param['bowtie2_type'])
    
    call.append('-N')
    call.append(param['bowtie2_N'])
    
    call.append('-D')
    call.append(param['bowtie2_D'])
    
    call.append('-R')
    call.append(param['bowtie2_R'])
    
    call.append('-L')
    call.append(param['bowtie2_L'])
    
    call.append('-i')
    call.append(param['bowtie2_i'])
    
    call.append('--rdg')
    call.append(param['bowtie2_rdg'])
    
    call.append('--rfg')
    call.append(param['bowtie2_rfg'])
    
    #index
    call.append('-x')
    call.append(param['bowtie2_index'])
    
    #number of processors
    call.append('-p')
    call.append(param['num_processors'])

    #if paired add second working file
    if param['paired']:
    	call.append('-1')
    	call.append(param['working_file'])
    	
    	call.append('-2')
        call.append(param['working_file2'])
    else:
    	call.append('-U')
    	call.append(param['working_file'])
    
    call.append('-S')
    call.append(outdir + 'accepted_hits.sam')

    param['file_handle'].write('CALL: '+' '.join(call)+'\n')
    output, error = subprocess.Popen(call, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    param['file_handle'].write(error)
    param['file_handle'].write(output)
    
    # Convert .sam to .bam and sort
    outfile = open(outdir + 'accepted_hits.bam','w')
    
    call2 = ['samtools']
    call2.append('view')
    call2.append('-bS')
    call2.append(outdir + 'accepted_hits.sam')
    
    param['file_handle'].write('CALL: '+' '.join(call2)+'\n')
    output, error = subprocess.Popen(call2, stdout = outfile, stderr=subprocess.PIPE).communicate()
    outfile.close()
    
    call3 = ['samtools']
    call3.append('sort')
    call3.append(outdir + 'accepted_hits.bam')
    call3.append(outdir + 'accepted_hits')
    
    param['file_handle'].write('CALL: '+' '.join(call3)+'\n')
    output, error = subprocess.Popen(call3, stdout = subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    param['file_handle'].write(error)
    param['file_handle'].write(output)
    
    # Remove SAM file
    call4 = ['rm']
    call4.append('-f')
    call4.append(outdir + 'accepted_hits.sam')
    
    param['file_handle'].write('CALL: '+' '.join(call4)+'\n')
    output, error = subprocess.Popen(call4, stdout = subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    param['file_handle'].write(error)
    param['file_handle'].write(output)

    #check if the run was successful
    if not os.path.exists(outdir+'accepted_hits.bam'):
        param['file_handle'].write('bowtie2 did not run successfully...')
        sys.exit(0)

    #wrap up and return the current workingfile
    MODULE_HELPER.wrapup_module(param, 
                                [outdir+'accepted_hits.bam'],
                                remove_intermediate=True)
    




