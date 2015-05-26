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

"""Wrapper to run the matched pairs script on all samples
"""
import hydra.module_helper
import os
import sys
import subprocess

def init(param):
    """Initializes all module specific parameter

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    hydra.module_helper.check_parameter(param,
                                                  key='match_pairs_exec',
                                                  dtype=str)


def run_match_pairs(param, infile, infile2, outfile, outfile2):
    """run the script that matches paired end sequencing reads

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    :Parameter infile: File that contains the first mates
    :Parameter infile2: File that contains the second mates
    :Parameter outfile: Where to write the cleaned up first mates
    :Parameter outfile2: Where to write the cleaned up second mates
    """

    call = [param['match_pairs_exec']]
    call.append(param[infile])
    call.append(param[infile2])
    call.append(outfile)
    call.append(outfile2)

    param['file_handle'].write(' '.join(call))
    output, error = subprocess.Popen(call,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE).communicate()

    param['file_handle'].write(error)
    param['file_handle'].write(output)

    with open(outfile+'.txt', 'w') as filehandle:
        filehandle.write(output)

    # Error handling
    if not os.path.exists(outfile2):
        param['file_handle'].write('ERROR: Could not find both output files')
        sys.exit(0)
    else:
        #check if the file integrity is alright
        call = ['gzip', '-t', outfile2] 
        output, error = subprocess.Popen(call,
                                         stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE).communicate()
        if 'unexpected end of file' in error:
            param['file_handle'].write('ERROR: File integrity corrupted, please rerun')
            sys.exit(0)
            
def main():
    """Main function that is run on each samples, which in turn calls the
    actual paired mate script that matches the mates
    """
    param = hydra.module_helper.initialize_module()

    #run match pairs
    outfile = (param['module_dir']+
               param['outstub']+
               '.clipped.matched.fastq.gz')
    outfile2 = (param['module_dir']+
                param['outstub']+
                '.clipped.matched.2.fastq.gz')

    run_match_pairs(param,
                    'working_file',
                    'working_file2',
                    outfile,
                    outfile2)
    hydra.module_helper.wrapup_module(param,
                                      [outfile, outfile2],
                                      remove_intermediate=True)

