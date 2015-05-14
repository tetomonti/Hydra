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

"""Star module
This module contains functions for initializing all tophat specific variables and
a wrapper that that runs tophat using those parameters on a single sample.
"""
import subprocess
import os
import sys
import hydra.module_helper
MODULE_HELPER = hydra.module_helper

def init(param):
    """Initialization function that checks the all relevant tophat parameters

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    MODULE_HELPER.check_parameter(param, key='star_exec', dtype=str)
    MODULE_HELPER.check_parameter(param, key='star_index', dtype=str, checkfile=True)
    MODULE_HELPER.check_parameter(param, key='outFilterType', dtype=str)
    MODULE_HELPER.check_parameter(param, key='outFilterMultimapNmax', dtype=str)
    MODULE_HELPER.check_parameter(param, key='alignSJoverhangMin', dtype=str)
    MODULE_HELPER.check_parameter(param, key='alignSJDBoverhangMin', dtype=str)
    MODULE_HELPER.check_parameter(param, key='outFilterMismatchNmax', dtype=str)
    MODULE_HELPER.check_parameter(param, key='outFilterMismatchNoverLmax', dtype=str)
    MODULE_HELPER.check_parameter(param, key='alignIntronMin', dtype=str)
    MODULE_HELPER.check_parameter(param, key='alignIntronMax', dtype=str)
    MODULE_HELPER.check_parameter(param, key='alignMatesGapMax', dtype=str)
    MODULE_HELPER.check_parameter(param,
                                  key='outputSAMtype',
                                  allowed=['BAM_SortedByCoordinate',
                                           'BAM_unsorted'],
                                  dtype=str)

def main():
    """Main function that is run on each samples, which in turn calls runs
    star on a sample.
    """

    param = MODULE_HELPER.initialize_module()
    #run create output directory
    outdir = param['module_dir']+param['outstub']+'/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    #build tophat call:
    call = [param['star_exec']]

    #add the directory where we built the star index
    call.append('--genomeDir')
    call.append(param['star_index'])

    #add the number of processors to use
    call.append('runThreadN')
    call.append(param['num_processors'])

    #add all the optional parameters
    call.append('--outFilterType')
    call.append(param['outFilterType'])
    call.append('--outFilterMultimapNmax')
    call.append(param['outFilterMultimapNmax'])
    call.append('--alignSJoverhangMin')
    call.append(param['alignSJoverhangMin'])
    call.append('--alignSJDBoverhangMin')
    call.append(param['alignSJDBoverhangMin'])
    call.append('--outFilterMismatchNmax')
    call.append(param['outFilterMismatchNmax'])
    call.append('--outFilterMismatchNoverLmax')
    call.append(param['outFilterMismatchNoverLmax'])
    call.append('--alignIntronMin')
    call.append(param['alignIntronMin'])
    call.append('--alignIntronMax')
    call.append(param['alignIntronMax'])
    call.append('--alignMatesGapMax')
    call.append(param['alignMatesGapMax'])
    if (param['outputSAMtype'] == 'BAM_SortedByCoordinate'):
        call.append('--outSAMtype')
        call.append('BAM')
        call.append('SortedByCoordinate')
        outfile = 'Aligned.sortedByCoord.out.bam'
    elif (param['outputSAMtype'] == 'BAM_unsorted'):
        call.append('--outSAMtype')
        call.append('BAM')
        call.append('Unsorted')
        outfile = 'Aligned.out.bam'
    else:
        outfile = 'Aligned.out.sam'
    #add the proper output directories
    call.append('--outFileNamePrefix')
    call.append(outdir)

    #specify whether the fastq files are zipped
    call.append('--readFilesCommand')
    if param['zipped_fastq']:
        call.append('gunzip')
        call.append('-c')
    else:
        call.append('UncompressionCommand')

    #adding the files we want to work on
    call.append('--readFilesIn')
    call.append(param['working_file'])

    #if paired add second working file
    if param['paired']:
        call.append(param['working_file2'])

    param['file_handle'].write('CALL: '+' '.join(call)+'\n')
    output, error = subprocess.Popen(call,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE).communicate()
    param['file_handle'].write(error)
    param['file_handle'].write(output)

    #check if the run was successful
    if not os.path.exists(outdir+'SJ.out.tab'):
        param['file_handle'].write('Star did not run successfully...')
        sys.exit(0)


    #wrap up and return the current workingfile
    MODULE_HELPER.wrapup_module(param, [outdir+outfile])




