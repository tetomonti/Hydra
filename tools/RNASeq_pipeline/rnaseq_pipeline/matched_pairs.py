"""Wrapper to run the matched pairs script on all samples
"""
import rnaseq_pipeline.module_helper
import os
import sys
import subprocess

def init(param):
    """Initializes all module specific parameter

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    rnaseq_pipeline.module_helper.checkParameter(param,
                                                 key='match_pairs_exec',
                                                 dType=str)


def run_match_pairs(param, infile, infile2, outfile, outfile2):
    """run the script that matches paired end sequencing reads

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    :Parameter infile: File that contains the first mates
    :Parameter infile2: File that contains the second mates
    :Parameter outfile: Where to write the cleaned up first mates
    :Parameter outfile2: Where to write the cleaned up second mates
    """

    call = param['match_pairs_exec']+' "'+param[infile]+ \
           '" "'+param[infile2]+'" '+outfile+' '+outfile2
    output, error = subprocess.Popen(call,
                                     shell=True,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE).communicate()

    param['file_handle'].write(error)
    param['file_handle'].write(output)

    with open(outfile+'.txt', 'w') as filehandle:
        filehandle.write(output)

    # Error handling
    if (not os.path.exists(outfile)) or os.stat(outfile).st_size < 1000:
        sys.exit(0)


def main():
    """Main function that is run on each samples, which in turn calls the \
    actual paired mate script that matches the mates
    """
    param = rnaseq_pipeline.module_helper.initialize_module()

    #run match pairs
    outfile = param['module_dir']+param['stub'][param['file_index']]+ \
              '.clipped.matched.fastq.gz'
    outfile2 = param['module_dir']+param['stub'][param['file_index']]+ \
               '.clipped.matched.2.fastq.gz'

    run_match_pairs(param,
                    'working_file',
                    'working_file2',
                    outfile,
                    outfile2)
    rnaseq_pipeline.module_helper.wrapup_module(param,
                                                [outfile, outfile2])

