"""Pipeline.py
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
import rnaseq_pipeline.helper
HELPER = rnaseq_pipeline.helper
import rnaseq_pipeline.fastqc
import rnaseq_pipeline.tophat
import rnaseq_pipeline.bamqc
import rnaseq_pipeline.cufflinks
import rnaseq_pipeline.matched_pairs
import rnaseq_pipeline.cutadapt
import rnaseq_pipeline.htseq
import rnaseq_pipeline.featureCount
import rnaseq_pipeline.star

import sys


def initialize_all(param):
    """this function calls the initialize functions of every module

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    rnaseq_pipeline.cutadapt.init(param)
    rnaseq_pipeline.matched_pairs.init(param)
    rnaseq_pipeline.fastqc.init(param)
    rnaseq_pipeline.tophat.init(param)
    rnaseq_pipeline.bamqc.init(param)
    rnaseq_pipeline.cufflinks.init(param)
    rnaseq_pipeline.htseq.init(param)
    rnaseq_pipeline.featureCount.init(param)
    rnaseq_pipeline.star.init(param)

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
            rnaseq_pipeline.cufflinks.finalize(param,
                                               input_files='count_files')

        if param['run_htseq']:
            HELPER.submit_job(param,
                              'run_htseq',
                              input_files='bam_files',
                              output_files='count_files')
            rnaseq_pipeline.htseq.finalize(param,
                                           input_files='count_files')

        if param['run_featureCount']:
            HELPER.submit_job(param,
                              'run_featureCount',
                              input_files='bam_files',
                              output_files='count_files')
            rnaseq_pipeline.featureCount.finalize(param,
                                                  input_files='count_files')

def report_all(param):
    """this function calls the reporting functions of every module

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    if param['aligner'] != 'skip':
        rnaseq_pipeline.fastqc.report(param,
                                      input_files='raw_files',
                                      header='FastQC results on the raw data')
        rnaseq_pipeline.fastqc.report(param,
                                      input_files='fastq_files',
                                      header='FastQC results after preprocessing')
    rnaseq_pipeline.bamqc.report(param)
    if param['run_cufflinks']:
        rnaseq_pipeline.cufflinks.report(param)
    if param['run_htseq']:
        rnaseq_pipeline.htseq.report(param)
    if param['run_featureCount']:
        rnaseq_pipeline.featureCount.report(param)



