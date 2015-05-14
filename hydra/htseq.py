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

"""HTSeq module
This module contains functions for initializing all HTSEQ specific variables and
a wrapper that that runs tophat using those parameters on a single sample.
In addition it also contains functions to extract and write statistics and
a wrapper that calls an R script
"""
import hydra.module_helper
MODULE_HELPER = hydra.module_helper
import hydra.helper
HELPER = hydra.helper
from hydra.r_scripts import get_script_path
import subprocess
import os

def init(param):
    """Initialization function that checks the all relevant tophat parameters

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    MODULE_HELPER.check_parameter(param, key='HTSeq_exec', dtype=str)
    MODULE_HELPER.check_parameter(param, key='sam_exec', dtype=str)
    MODULE_HELPER.check_parameter(param, key='HTSeq_s', dtype=str)
    MODULE_HELPER.check_parameter(param, key='HTSeq_t', dtype=str)
    MODULE_HELPER.check_parameter(param,
                                  key='HTSeq_m',
                                  dtype=str,
                                  allowed=['union',
                                           'intersection-strict',
                                           'intersection-nonempty'])
    MODULE_HELPER.check_parameter(param, key='HTSeq_id', dtype=str)
    MODULE_HELPER.check_parameter(param, key='HTSeq_gft', dtype=str, checkfile=True)
    MODULE_HELPER.check_parameter(param, key='Rscript_exec', dtype=str)

def create_eset(count_file, pheno_file, param):
    """Wrapper that calls an R script which creates a Bioconductor ExpressionSet

    :Parameter count_file: location of the count file
    :Parameter pheno_file: location of the phenotype file
    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    HELPER.writeLog('Creating ESet ... \n', param)
    #create a Bioconductor ExpresionSet
    call = [param['Rscript_exec']]
    call.append(get_script_path('createRawCountESet.R'))
    call = call + ['-c', count_file]
    call = call + ['-a', pheno_file]
    call = call + ['-o', param['working_dir']]
    call = call + ['-s', 'htseq']
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
    param['report'].write('<br><a href="htseq/htseq_pca.html">PCA on normalized samples</a>')

def process_stat_files(param):
    """Copies all relevant files into the report directory and also extracts
    the total number of reads from the bamqc output

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    #if there is no htseq directory in the report make one
    htseq_dir = param['working_dir']+'report/htseq/'
    if not os.path.exists(htseq_dir):
        os.makedirs(htseq_dir)

    #get the files that are actually in the output directory
    call = ['cp', '-R']
    call.append(param['working_dir']+'results/htseq/')
    call.append(param['working_dir']+'report/')
    _, _ = subprocess.Popen(call,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE).communicate()

    htseq_file = param['working_dir']+'results/htseq/htseq_stats.txt'
    #extract table
    table = []
    filehandle = open(htseq_file)
    #header
    table.append(filehandle.readlines()[0].rstrip().split('\t'))
    table[0] = table[0][1:]
    filehandle.close()

    #total number of aligned reads
    tot_reads = ['Total number of aligned reads']+param['bam_qc']['total_aligned_reads']
    table.append(tot_reads)

    filehandle = open(htseq_file)
    for line in filehandle.readlines()[1:]:
        cur_line = line.rstrip().split('\t')
        perc = ([cur_line[0]]+
                MODULE_HELPER.get_percentage(cur_line[1:],
                                             tot_reads[1:],
                                             len(cur_line)-1))
        table.append(perc)

    filehandle.close()
    return table



def report(param):
    """Function that writes all HTSeq related statistics into the html report

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """

    out_file = param['working_dir']+'deliverables/htseq_raw_counts.txt'
    #report only if there were actually results
    if os.path.exists(out_file):
        param['report'].write('<center><br><br><br><br><h2>HTSeq statistics</h2>')
        table = process_stat_files(param)
        MODULE_HELPER.write_html_table(param,
                                       table,
                                       out=param['report'],
                                       cell_width=80,
                                       fcol_width=150,
                                       deg=315)
        param['report'].write('<a href="htseq/htseq_stats.txt">'+
                              'HTSeq statistics as tab delimited txt file</a>')
        #create an eSet:
        create_eset(out_file,
                    param['pheno_file'],
                    param)



def finalize(param, input_files='count_files'):
    """This function is run after HTSeq is run on each sample. It collects all results
    and puts them into a file

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    :Parameter input_files: flag that indicates the input files
    """
    HELPER.writeLog('Collecting HTSeq raw counts ... \n', param)
    #extracts the counts from the htseq output
    import csv

    #check which of these files are actually available
    working_files = [iFile for iFile in param[input_files] if iFile != '']

    if len(working_files) > 0:
        #get gene annotation
        csv_file = open(working_files[0])
        csv_reader = csv.reader(csv_file, delimiter='\t')
        counts = [row[0] for row in csv_reader]
        csv_file.close()

        #get all the expression values
        header = 'ID'
        for idx in range(param['num_samples']):
            if param[input_files] != '':
                header = header+'\t'+param['stub'][idx]
                csv_file = open(param[input_files][idx])
                csv_reader = csv.reader(csv_file, delimiter='\t')
                i = 0
                for row in csv_reader:
                    counts[i] = counts[i]+'\t'+row[1]
                    i += 1
                csv_file.close()

        #output the file
        out_file = param['working_dir']+'deliverables/htseq_raw_counts.txt'
        out_handle = open(out_file, 'w')
        out_handle.write(header+'\n')
        #drop the last 5 lines since these provide only a summary statistic
        for i in range(len(counts)-5):
            out_handle.write(counts[i]+'\n')
        out_handle.close()

        #output_phenotype_file
        HELPER.writeLog('Writing phenotype data ... \n', param)
        MODULE_HELPER.output_sample_info(param)

        #write summary stats
        out_handle = open(param['working_dir']+
                          'results/htseq/htseq_stats.txt',
                          'w')
        out_handle.write(header+'\n')
        for i in range(len(counts)-5, len(counts)):
            out_handle.write(counts[i]+'\n')
        out_handle.close()

    else:
        print 'HTseq was not run successfully on any of the files..\n'



def main():
    """Main function that is run on each samples, which in turn calls runs
    htseq on a sample.
    """
    import sys
    param = MODULE_HELPER.initialize_module()

    #build htseq-count call:
    call1 = [param['sam_exec'], 'view', param['working_file']]
    call2 = [param['HTSeq_exec']]
    call2 = call2 + ['-s', param['HTSeq_s']]
    call2 = call2 + ['-t', param['HTSeq_t']]
    call2 = call2 + ['-i', param['HTSeq_id']]
    call2 = call2 + ['-m', param['HTSeq_m']]
    call2 = call2 + ['-', param['HTSeq_gft']]

    #function calls
    param['file_handle'].write('Pipe CALL 1: '+' '.join(call1)+'\n')
    param['file_handle'].write('Pipe CALL 2: '+' '.join(call2)+'\n')

    process1 = subprocess.Popen(call1, stdout=subprocess.PIPE)
    process2 = subprocess.Popen(call2, stdin=process1.stdout, stdout=subprocess.PIPE)
    process1.stdout.close()
    output, error = process2.communicate()

    #error handling
    if output == '':
        param['file_handle'].write(error+'\n')
        sys.exit(0)

    #write output
    outfile = (param['module_dir']+
               param['outstub']+
               '.txt')
    handle = open(outfile, 'w')
    handle.write(output)
    handle.close()

    #wrap up and return the current workingfile
    MODULE_HELPER.wrapup_module(param, [outfile])



