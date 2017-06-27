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
import hydra_pkg.helper
HELPER = hydra_pkg.helper
import hydra_pkg.fastqc
import hydra_pkg.tophat
import hydra_pkg.bamqc
import hydra_pkg.cufflinks
import hydra_pkg.matched_pairs
import hydra_pkg.cutadapt
import hydra_pkg.htseq
import hydra_pkg.featureCount
import hydra_pkg.star
import hydra_pkg.bowtie2

import sys


def initialize_all(param):
    """this function calls the initialize functions of every module

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    hydra_pkg.cutadapt.init(param)
    hydra_pkg.fastqc.init(param)
    hydra_pkg.bamqc.init(param)
    hydra_pkg.cufflinks.init(param)
    hydra_pkg.htseq.init(param)
    hydra_pkg.featureCount.init(param)

    if param['aligner'] == 'tophat':
        hydra_pkg.tophat.init(param)
    elif param['aligner'] == 'star':    
        hydra_pkg.star.init(param)
    elif param['aligner'] == 'bowtie2':
        hydra_pkg.bowtie2.init(param)

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

        if not param['skip_trimming']:
            HELPER.submit_job(param,
                              'run_cutadapt',
                              input_files='raw_files',
                              output_files='fastq_files')
            HELPER.submit_job(param,
                              'run_fastqc',
                              input_files='fastq_files')

            #if indicated remove samples that failed the QC
            if param['remove_failed']:
                hydra_pkg.fastqc.remove_failed(param, 
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
            elif param['aligner'] == 'star':
                #running the aligner
                HELPER.submit_job(param,
                                  'run_star',
                                  input_files='fastq_files',
                                  output_files='bam_files',
                                  cores=param['qsub_num_processors'],
                                  mem_free='32G')
            elif param['aligner'] == 'bowtie2':
                #running the aligner
                HELPER.submit_job(param,
                                  'run_bowtie2',
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
            hydra_pkg.cufflinks.finalize(param,
                                         input_files='count_files')

        if param['run_htseq']:
            HELPER.submit_job(param,
                              'run_htseq',
                              input_files='bam_files',
                              output_files='count_files')
            hydra_pkg.htseq.finalize(param,
                                     input_files='count_files')

        if param['run_featureCount']:
            HELPER.submit_job(param,
                              'run_featureCount',
                              input_files='bam_files',
                              output_files='count_files')
            hydra_pkg.featureCount.finalize(param,
                                            input_files='count_files')

def report_all(param):
    """this function calls the reporting functions of every module

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    if param['aligner'] != 'skip':
        hydra_pkg.fastqc.report(param,
                                input_files='raw_files',
                                header='FastQC results on the raw data')
        if not param['skip_trimming']:
            hydra_pkg.fastqc.report(param,
                                    input_files='fastq_files',
                                    header='FastQC results after preprocessing')
    hydra_pkg.bamqc.report(param)
    if param['run_cufflinks']:
        hydra_pkg.cufflinks.report(param)
    if param['run_htseq']:
        hydra_pkg.htseq.report(param)
    if param['run_featureCount']:
        hydra_pkg.featureCount.report(param)



