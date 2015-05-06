"""Cutadapt module
This module contains functions for running the cutadapt tool
"""

import rnaseq_pipeline.module_helper
MODULE_HELPER = rnaseq_pipeline.module_helper
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


def run_cutadapt(param, infile, outfile, adapter):
    """Runs cutadapt on a file that remove an adapter sequence

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    :Parameter infile: input filename
    :Parameter outfile: output filename
    :Parameter adapter: adapter sequence

    """
    call = [param['cutadapt_exec']]
    call = call + ['-m', param['cutadapt_m']]
    call = call + ['-a', adapter]
    call = call + ['-o', outfile]
    call.append(param[infile])

    param['file_handle'].write(' '.join(call))
    output, error = subprocess.Popen(call,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE).communicate()

    param['file_handle'].write(error)
    param['file_handle'].write(output)

    with open(outfile+'.txt', 'w') as filehandle:
        filehandle.write(output)

    # Error handling
    if (not os.path.exists(outfile)) or os.stat(outfile).st_size < 10000:
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
    run_cutadapt(param,
                 'working_file',
                 outfile,
                 param['cutadapt_first_adapter'])

    if not param['paired']:
        MODULE_HELPER.wrapup_module(param, [outfile])
    #calling it on the second fastq file if it is paired
    else:
        outfile2 = (param['module_dir']+
                    param['outstub']+
                    '.clipped.2.fastq.gz')
        run_cutadapt(param,
                     'working_file2',
                     outfile2,
                     param['cutadapt_second_adapter'])
        MODULE_HELPER.wrapup_module(param, [outfile, outfile2])
