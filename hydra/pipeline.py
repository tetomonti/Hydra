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

"""pipeline.py
This script specifies the flow of the pipeline and indicates in which order
the single modules are strung together. Since all modules are programmed fairly
dumb, i.e. they do not know of modules that come before or after it enables
a simple way of changing the flow of the pipeline and  is also the only place
other than the parameter file, which has to be modified in order to include
a new module, e.g. another QC or another aligner.

In addition to the workflow, all module specific parameters are initialized here
before the pipeline starts running and the all module specific reporting functions
are executed as well.
"""
import hydra.helper
HELPER = hydra.helper
import hydra.fastqc
import hydra.tophat
import hydra.bamqc
import hydra.cufflinks
import hydra.matched_pairs
import hydra.cutadapt
import hydra.htseq
import hydra.featureCount
import hydra.star

import sys


def initialize_all(param):
    """this function calls the initialize functions of every module

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    hydra.cutadapt.init(param)
    hydra.matched_pairs.init(param)
    hydra.fastqc.init(param)
    hydra.tophat.init(param)
    hydra.bamqc.init(param)
    hydra.cufflinks.init(param)
    hydra.htseq.init(param)
    hydra.featureCount.init(param)
    hydra.star.init(param)

def run_all(param):
    """this function defines the workflow of the pipeline

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    if param['aligner'] == 'skip':
        param['bam_files'] = param['raw_files'][:]
    else:
        #preprocessing fastq file
        HELPER.submit_job(param, 'run_fastqc', input_files='raw_files')
        param['fastq_files'] = param['raw_files'][:]
        if param['paired']:
            param['fastq_files2'] = param['raw_files2'][:]
        HELPER.submit_job(param,
                          'run_cutadapt',
                          input_files='raw_files',
                          output_files='fastq_files')
        if param['paired']:
            HELPER.submit_job(param,
                              'run_matched_pairs',
                              input_files='fastq_files',
                              output_files='fastq_files')
        HELPER.submit_job(param,
                          'run_fastqc',
                          input_files='fastq_files')

        #do alignment if it's not just a fastqc run
        if not param['QC_and_trim_only']:
            if param['aligner'] == 'tophat':
                #running the aligner
                HELPER.submit_job(param,
                                  'run_tophat',
                                  input_files='fastq_files',
                                  output_files='bam_files',
                                  cores=param['qsub_num_processors'])
            if param['aligner'] == 'star':
                #running the aligner
                HELPER.submit_job(param,
                                  'run_star',
                                  input_files='fastq_files',
                                  output_files='bam_files',
                                  cores=param['qsub_num_processors'],
                                  mem_free='32G')
            else:
                HELPER.writeLog('The selected aligner does not exist.', param)
                sys.exit(0)

    if not param['QC_and_trim_only']:
        #Bamqc
        HELPER.submit_job(param,
                          'run_bamqc',
                          input_files='bam_files')

        #Getting the counts:
        if param['run_cufflinks']:
            HELPER.submit_job(param,
                              'run_cufflinks',
                              input_files='bam_files',
                              output_files='count_files')
            hydra.cufflinks.finalize(param,
                                               input_files='count_files')

        if param['run_htseq']:
            HELPER.submit_job(param,
                              'run_htseq',
                              input_files='bam_files',
                              output_files='count_files')
            hydra.htseq.finalize(param,
                                           input_files='count_files')

        if param['run_featureCount']:
            HELPER.submit_job(param,
                              'run_featureCount',
                              input_files='bam_files',
                              output_files='count_files')
            hydra.featureCount.finalize(param,
                                                  input_files='count_files')

def report_all(param):
    """this function calls the reporting functions of every module

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    if param['aligner'] != 'skip':
        hydra.fastqc.report(param,
                                      input_files='raw_files',
                                      header='FastQC results on the raw data')
        hydra.fastqc.report(param,
                                      input_files='fastq_files',
                                      header='FastQC results after preprocessing')
    hydra.bamqc.report(param)
    if param['run_cufflinks']:
        hydra.cufflinks.report(param)
    if param['run_htseq']:
        hydra.htseq.report(param)
    if param['run_featureCount']:
        hydra.featureCount.report(param)



