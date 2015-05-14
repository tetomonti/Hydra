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

"""Tophat module
This module contains functions for initializing all tophat specific variables and
a wrapper that that runs tophat using those parameters on a single sample.
"""

import hydra.module_helper
MODULE_HELPER = hydra.module_helper

def init(param):
    """Initialization function that checks the all relevant tophat parameters

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    MODULE_HELPER.check_parameter(param, key='tophat_exec', dtype=str)
    MODULE_HELPER.check_parameter(param, key='tophat_index', dtype=str)
    MODULE_HELPER.check_parameter(param, key='tophat_qual', dtype=str, optional=True)
    MODULE_HELPER.check_parameter(param, key='tophat_N', dtype=str)
    MODULE_HELPER.check_parameter(param, key='tophat_gap_length', dtype=str)
    MODULE_HELPER.check_parameter(param, key='tophat_edit_dist', dtype=str)
    MODULE_HELPER.check_parameter(param, key='mate_inner_dist', dtype=str)
    MODULE_HELPER.check_parameter(param, key='mate_std_dev', dtype=str)

def main():
    """Main function that is run on each samples, which in turn calls runs
    tophat on a sample.
    """
    import subprocess
    import sys
    import os
    param = MODULE_HELPER.initialize_module()
    #run create output directory
    outdir = param['module_dir']+param['outstub']+'/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    #build tophat call:
    call = [param['tophat_exec']]

    #resume option
    if param['resume_module']:
        call.append('-R')
        call.append(outdir)

    else:
        #optional quality value
        if param['tophat_qual'] != 'none':
            call.append(param['tophat_qual'])

        #add paired parameters if paired
        if param['paired']:
            call.append('--mate-inner-dist')
            call.append(param['mate_inner_dist'])
            call.append('--mate-std-dev')
            call.append(param['mate_std_dev'])

        call.append('-o')
        call.append(outdir)
        call.append('-p')
        call.append(param['num_processors'])
        call.append('-N')
        call.append(param['tophat_N'])
        call.append('--read-gap-length')
        call.append(param['tophat_gap_length'])
        call.append('--read-edit-dist')
        call.append(param['tophat_edit_dist'])
        call.append(param['tophat_index'])
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

    #error handling - checking if tophat logfile exists and if it does,
    #check if tophat was run successful
    tophat_logfile = outdir + 'logs/tophat.log'

    if not os.path.exists(tophat_logfile):
        sys.exit(0)
    else:
        log = open(tophat_logfile)

        lines_end = []
        for line in log.readlines():
            if ' Run complete: ' in line.rstrip():
                lines_end.append(line)
        log.close()
        if len(lines_end) == 0:
            sys.exit(0)

    #wrap up and return the current workingfile
    MODULE_HELPER.wrapup_module(param, [outdir+'accepted_hits.bam'])




