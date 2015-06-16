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

"""Cutadapt module
This module contains functions for running the cutadapt tool
"""

import hydra_pkg.module_helper as MODULE_HELPER
import os
import sys
import subprocess

def init(param):
    """Initialization function, that checks if the bamqc_script that is run
    on every single samples is available

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    MODULE_HELPER.check_parameter(param, key='cutadapt_exec', dtype=str)
    MODULE_HELPER.check_parameter(param, key='cutadapt_first_adapter', dtype=str)
    if param['paired']:
        MODULE_HELPER.check_parameter(param, key='cutadapt_second_adapter', dtype=str)
    MODULE_HELPER.check_parameter(param, key='cutadapt_m', dtype=str)
    MODULE_HELPER.check_parameter(param, key='cutadapt_q_end', dtype=str)
    MODULE_HELPER.check_parameter(param, key='cutadapt_q_start', dtype=str)
    MODULE_HELPER.check_parameter(param, key='cutadapt_quality', dtype=str)


def run_cutadapt(param, outfile, outfile2=''):
    """Runs cutadapt on a file that remove an adapter sequence

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    :Parameter infile: input filename
    :Parameter outfile: output filename
    """
    call = [param['cutadapt_exec']]
    call = call + ['-m', param['cutadapt_m']]
    call = call + ['-q', param['cutadapt_q_start']+','+param['cutadapt_q_end']]
    call = call + ['-a', param['cutadapt_first_adapter']]
    if param['paired']:
        call = call + ['-A', param['cutadapt_second_adapter']]    
    call = call + ['-o', outfile]
    if param['paired']:
        call = call + ['-p', outfile2]  
    call.append(param['working_file'])
    if param['paired']:
        call.append(param['working_file2'])

    param['file_handle'].write(' '.join(call))
    output, error = subprocess.Popen(call,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE).communicate()

    param['file_handle'].write(error)
    param['file_handle'].write(output)

    with open(outfile+'.txt', 'w') as filehandle:
        filehandle.write(output)


    # Error handling
    if not os.path.exists(outfile):
        param['file_handle'].write('ERROR: Could not find the output file')
        sys.exit(0)
    else:
        #check if the file integrity is alright
        call = ['gzip', '-t', outfile] 
        output, error = subprocess.Popen(call,
                                         stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE).communicate()
        if 'unexpected end of file' in error:
            param['file_handle'].write('ERROR: File integrity corrupted, please rerun')
            sys.exit(0)
            

def main():
    """Main function that is run on each samples, which in turn calls the
    cutadapt tool
    """
    param = MODULE_HELPER.initialize_module()

    #run fastqc
    outfile = (param['module_dir']+
               param['outstub']+
               '.clipped.fastq.gz')

    if not param['paired']:
        run_cutadapt(param, outfile)
        MODULE_HELPER.wrapup_module(param, [outfile])
    else:
        outfile2 = (param['module_dir']+
                    param['outstub']+
                    '.clipped.2.fastq.gz')
        run_cutadapt(param, outfile, outfile2)
        MODULE_HELPER.wrapup_module(param, [outfile, outfile2])
