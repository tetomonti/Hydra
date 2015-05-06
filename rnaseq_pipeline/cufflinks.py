"""Cufflinks module
This module contains functions for initializing all cufflinks specific variables and
a wrapper that that runs cufflinks using those parameters on a single sample,
as well as functions for reporting the results
"""

import rnaseq_pipeline.module_helper
MODULE_HELPER = rnaseq_pipeline.module_helper
import rnaseq_pipeline.helper
HELPER = rnaseq_pipeline.helper
from rnaseq_pipeline.r_scripts import get_script_path
import os, subprocess

def init(param):
    """Initialization function that checks the all relevant tophat parameters
    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    MODULE_HELPER.check_parameter(param, key='cufflinks_exec', dtype=str)
    MODULE_HELPER.check_parameter(param, key='cufflinks_G', dtype=str, checkfile=True)
    MODULE_HELPER.check_parameter(param, key='cufflinks_b', dtype=str, checkfile=True)
    MODULE_HELPER.check_parameter(param, key='cufflinks_library_type', dtype=str)
    MODULE_HELPER.check_parameter(param, key='cufflinks_compatible_hits', dtype=str)
    MODULE_HELPER.check_parameter(param, key='cufflinks_total_hits', dtype=str)
    MODULE_HELPER.check_parameter(param, key='cufflinks_N', dtype=str)
    MODULE_HELPER.check_parameter(param, key='cufflinks_u', dtype=str)

def report(param):
    """This function runs PCA and writes the corresponding link in the html report

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """

    param['report'].write('<center><br><br><h2>Cufflinks results</h2>')
    #make cufflinks directory
    cufflinks_dir = param['working_dir']+'report/cufflinks/'
    if not os.path.exists(cufflinks_dir):
        os.makedirs(cufflinks_dir)

    #output_phenotype_file
    MODULE_HELPER.output_sample_info(param)

    #run R script that creates a PCA
    HELPER.writeLog('Running PCA ... \n', param)
    counts = param['working_dir']+'deliverables/cufflinks_counts_fpkm.txt'
    call = [param['Rscript_exec']]
    call.append(get_script_path('cufflinks_pca.R'))
    call = call + ['-c', counts]
    call = call + ['-a', param['pheno_file']]
    call = call + ['-o', param['working_dir']]
    call = call + ['-s', 'cufflinks']
    if param['paired']:    
        paired = 'TRUE'
    else:
        paired = 'FALSE'
    call = call + ['-p', paired]

    output, error = subprocess.Popen(call,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE).communicate()
    HELPER.writeLog(output, param)
    HELPER.writeLog(error, param)
    param['report'].write('<a href="cufflinks/cufflinks_pca.html">PCA</a>')



def finalize(param, input_files='count_files'):
    """
    This function collects Cufflinks results and generates a text file with the
    gene identifiers and their counts

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    :param input_files: working files produced after running Cufflinks

    """
    HELPER.writeLog('Collecting Cufflinks FPKM data ... \n', param)
    #extracts the counts from the cufflinks output
    import csv
    #check which of these files are actually available
    working_files = [iFile for iFile in param[input_files] if iFile != '']
    if len(working_files) > 0:
        #get gene annotation
        csv_file = open(working_files[0])
        csv_reader = csv.reader(csv_file, delimiter='\t')
        counts = [row[0]+'\t'+row[4] for row in csv_reader]
        csv_file.close()
        #get all the expression values
        header = 'ENS_ID\tSymbol'
        for idx in range(param['num_samples']):
            if param[input_files] != '':
                header = header + '\t' + param['stub'][idx]
                csv_file = open(param[input_files][idx])
                csv_reader = csv.reader(csv_file, delimiter='\t')
                i = 0
                for row in csv_reader:
                    counts[i] = counts[i] + '\t' + row[9]
                    i += 1
                csv_file.close()
        #output the file
        out_handle = open(param['working_dir']+'deliverables/cufflinks_counts_fpkm.txt', 'w')
        out_handle.write(header+'\n')
        for i in range(len(counts)):
            out_handle.write(counts[i]+'\n')
        out_handle.close()
    else:
        print 'Cufflinks was not run successfully on any of the files..\n'

def main():
    """ Main function that runs Cufflinks call
    """
    import sys
    param = MODULE_HELPER.initialize_module()

    #run create output directory
    outdir = param['module_dir']+param['outstub']+'/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    #build cufflinks call:
    call = [param['cufflinks_exec']]
    #optional parameters
    if param['cufflinks_compatible_hits'] == 'active':
        call.append('--compatible-hits-norm')
    if param['cufflinks_N'] == 'active':
        call.append('-N')
    if param['cufflinks_u'] == 'active':
        call.append('-u')
    if param['cufflinks_total_hits'] == 'active':
        call.append('--total-hits-norm')

    #finish call
    call.append('--no-update-check')
    call.append('--library-type')
    call.append(param['cufflinks_library_type'])
    call = call + ['-o' + outdir]
    call = call + ['-p' + param['num_processors']]
    call = call + ['-G' + param['cufflinks_G']]
    call.append(param['working_file'])

    param['file_handle'].write('CALL: '+' '.join(call)+'\n')
    output, _ = subprocess.Popen(call,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE).communicate()
    param['file_handle'].write(output)

    #printing the error blows the log file up to several megabytes
    #param['file_handle'].write(error)
    #error handling
    if not os.path.exists(outdir+'genes.fpkm_tracking'):
        param['file_handle'].write('Error there was no cuffinks output\n')
        sys.exit(0)
    #wrap up and return the current workingfile
    MODULE_HELPER.wrapup_module(param, [outdir+'genes.fpkm_tracking'])