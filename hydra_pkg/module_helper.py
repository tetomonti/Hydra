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

"""Module helper
Contains all helper functions that are relevant for the wrappers that run
an external tool on each samples.
"""

import hydra_pkg.helper as HELPER
import os
import sys
import matplotlib.pyplot as plt
import subprocess
import getopt
import json
import re
from hydra_pkg.r_scripts import get_script_path

def check_parameter(param, key, dtype, allowed=[], checkfile=False, optional=False):
    """generic function that checks if a parameter was in the parameter file,
    casts to the right data type and if the parameter is a file/ directory
    checks if it actually exists

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    :Parameter key: key for the entry in the dictionary that we want to check
    :Parameter dType: data type that the value is cast to
    :Parameter allowed: values that are allowed as input value
    :Parameter checkFile: indicates that the value is a file or directory
                          which should be checked
    :Parameter optional: indicates a parameter is optional, if it doesn't exits
                         it will be created
    """
    HELPER.check_parameter(param=param,
                           key=key,
                           dtype=dtype,
                           allowed=allowed,
                           checkfile=checkfile,
                           optional=optional)


def get_percentage(number1, number2, ntotal):
    """generic function that checks if a parameter was in the parameter file,
    casts to the right data type and if the parameter is a file/directory
    checks if it actually exists

    :Parameter number1: array of numerators
    :Parameter number2: array of denominators
    :Parameter ntotal: length of the arrays
    :rtype: array with percentages as strings with ntotal length
    """
    if type(number1) is float:
        number1 = [number1]
    if type(number2) is float:
        number2 = [number2]
    percent = [0.0]*ntotal
    
    for idx in range(ntotal):
        percent[idx] = round(float(number1[idx])/float(number2[idx])*100, 2)
    return [str(round(pc, 1))+'%' for pc in percent]

def divide(num1, num2, ntotal):
    """helper function that returns percentages from 2 arrays

    :Parameter number1: array of numerators
    :Parameter number2: array of denominators
    :Parameter ntotal: length of the arrays
    :rtype: array with percentages as strings with ntotal length
    """    
    percent = [0.0]*ntotal    
    for idx in range(ntotal):
        percent[idx] = float(num1[idx])/float(num2[idx])*100
    return percent

def initialize_module():
    """Generic function that is executed by every module before it runs on a
    sample. Provides information on which file to work with, number of cores,
    current working directory, handles the log file opening and checks if
    the pipeline is run in resume mode, so it skips the entire module run
    if this step has been completed successfully.
    """

    working_dir = './'
    #Check arguments
    if len(sys.argv) < 5:
        print 'ERROR: Specify the index of the file the parameter should be run on.'
        sys.exit(0)
    optlist, _ = getopt.getopt(sys.argv[1:], 'i:n:d:')
    for opt in optlist:
        if opt[0] == '-i':
            file_index = opt[1]
        if opt[0] == '-n':
            num_processors = opt[1]
        if opt[0] == '-d':
            working_dir = opt[1]
    print '\n'
    print '###########################'
    print sys.argv
    print working_dir


    #Read and initialize parameters
    with open(working_dir+'results/parameters.json') as filehandle:
        param = json.load(filehandle)
    param['file_index'] = int(file_index)
    param['num_processors'] = num_processors
    param['outstub'] = param['stub'][param['file_index']]

    #use the input files that were specified in the pipeline call
    param['working_file'] = param[param['input_files']][param['file_index']]
    if param['paired'] and param['input_files']+'2' in param:
        param['working_file2'] = param[param['input_files']+'2'][param['file_index']]

    #name of the log file
    log_file = (param['working_dir']+
                'results/log/'+
                param['stub'][param['file_index']]+
                '.log')

    #check if it is not a clean run and if module already completed
    param['resume_module'] = False
    if not param['clean_run']:
        param['file_handle'] = open(log_file)
        #check if the module already finished
        lines_end = False
        for line in param['file_handle'].readlines():
            if 'ENDING %s |' %(param['current_flag']) in line.rstrip():
                lines_end = True
        if lines_end:
            param['file_handle'].close()
            #open file for writing
            param['file_handle'] = open(log_file, 'a')
            param['file_handle'].write(param['current_flag']+
                                       ' module already run on this file .. SKIPPING\n')
            param['file_handle'].close()
            sys.exit(0)
        #check if the module was started, but not finished the enable resuming

        for line in param['file_handle'].readlines():
            if 'STARTING %s |' %(param['current_flag']) in line.rstrip():
                param['resume_module'] = True

    #start process log
    param['file_handle'] = open(log_file, 'a')
    param['file_handle'].write('STARTING '+param['current_flag']+'\n')
    return param

def output_phenotype(param, pheno_file):
    """Writes out a phenotype file

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    :Parameter pheno_file: filename of the output file
    """

    #write sample info only for the files that successfully completed
    header = param['raw_file_header']
    index = 0
    out = open(pheno_file, "w")
    filehandle = open(param['raw_filenames'], 'r')
    for line in filehandle:
        #take header into account
        if header:
            out.write(line)
            header = False
        else:
            #write only samples that successfully finished
            if param['run_log'][-1][index]:
                out.write(line)
            index += 1
    out.close()
    filehandle.close()

def output_sample_info(param):
    """Writes a phenotype file of the samples that actually made it through HTSeq
    into the report directory

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    param['pheno_file'] = param['working_dir']+'deliverables/sample_info.txt'
    output_phenotype(param, param['pheno_file'])


def is_in_raw_files(raw_files, current_file):
    """Function to check whether the fiel we want to delete is a raw file

    :Parameter raw_files: list of raw file loaction
    :Parameter current_file: file that we want to delete
    """
    is_in_raw = False
    for temp in raw_files:
        if temp == current_file:
            is_in_raw = True
    return is_in_raw
    
def wrapup_module(param, new_working_file=[], remove_intermediate=False):
    """Function to wrap up a module run. Writes the ending flag which is used to
    identify if the module was completed correctly. Also closes the log file handle.
    And finally also sets the working file pointer to the output of the module
    if the module has output files that are used in the next step.

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    :Parameter new_working_file: pointer towards a potential new working file
    """
 
    #remove intermediate files if specified
    if (remove_intermediate):
        #this is just a defensive check to make sure the raw files 
        #cannot be touched even if the flag is specified wrongly and someone 
        #changes the order the modules are run        
        if not is_in_raw_files(param['raw_files'], param['working_file']):
            os.remove(param['working_file'])
        else:
            param['file_handle'].write('WARNING: A module tried to delete a'+
                                       ' raw file. SKIPPING the removal, please'+
                                       ' fix the source code!')
        if param['paired']:
            if not is_in_raw_files(param['raw_files'], param['working_file']):
                os.remove(param['working_file2'])

    #end process log
    param['file_handle'].write('ENDING '+param['current_flag']+' | ')

   #if there was an actual output file specified
    if len(new_working_file) > 0:
        param['file_handle'].write(';'.join([w for w in new_working_file]))
    param['file_handle'].write('\n\n')
    param['file_handle'].close()
    
    
def plot_count_overview(param, stub):
    """Function that plot the overview boxplots for featureCounts and HTSeq

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    :Parameter stub: stub that indicates which module we are working with
    """
    
    #extract table
    overview = []
    labels = []
 
    #total number of aligned reads:
    tot_reads = param['bam_qc']['total_aligned_reads']

    #if we have paired end reads we need to divide by 2 to get fragments rather than reads
    if param['paired']:
        tot_reads = [reads / 2 for reads in tot_reads]

    #read in stats file
    stats_file = param['working_dir']+'results/'+stub+'/'+stub+'_stats.txt'
    filehandle = open(stats_file)
    for line in filehandle.readlines()[1:]:
        cur_line = line.rstrip().split('\t')
        labels.append(re.sub(r'_',' ',cur_line[0]))
        perc = (divide(cur_line[1:],
                tot_reads,
                len(cur_line)-1))
        overview.append(perc)
    filehandle.close()

    #make the first plot out of the first 2:
    fig, ax = plt.subplots()
    fig.set_size_inches(9, len(labels) / 2 + 0.5)
    bp = ax.boxplot(overview, patch_artist=True, vert=False)

    #change coloring
    for box in bp['boxes']:
        box.set( color='#7570b3', linewidth=2)
        box.set( facecolor = '#999999' )

    #change caps
    for cap in bp['caps']:
        cap.set(color='#7570b3', linewidth=2)

    #change outliers
    for flier in bp['fliers']:
        flier.set(marker='o', color='#ff0000', alpha=0.5)

    ax.set_yticklabels(labels) 
    ax.set_xlim(-5,105)

    #put it into the report
    filename = param['working_dir']+'report/'+stub+'/overview.png'
    fig.savefig(filename, bbox_inches='tight')
    param['report'].write('<img src="'+stub+'/overview.png" ' +
                          'alt="overview"><br><br>\n')




def create_eset(count_file, pheno_file, param, stub):
    """Wrapper that calls an R script which creates a Bioconductor ExpressionSet

    :Parameter count_file: location of the count file
    :Parameter pheno_file: location of the phenotype file
    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    :Parameter stub: stub indicating the module we are working with    
    """
    HELPER.writeLog('Creating ESet ... \n', param)
    #create a Bioconductor ExpresionSet
    call = [param['Rscript_exec']]
    call.append(get_script_path('createRawCountESet.R'))
    call = call + ['-c', count_file]
    call = call + ['-a', pheno_file]
    call = call + ['-o', param['working_dir']]
    call = call + ['-s', stub]
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
    param['module_report'].write('<br><a href="' + stub + '_pca.html">PCA on normalized samples</a>')
    param['module_report'].write('<br><center><h3>Boxplot of counts in log2 space</h3>')
    param['module_report'].write('<img src="' + stub + '_boxplot.png"' +
                                ' alt="Boxplot of ' + stub +' counts"><br><br>\n')
                                
                                
def create_sub_report(param, out_file, table, stub, title):
    """Separate report for all the htseq/featureCount results in detail

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    :Parameter out_file: report html filehandle
    :Parameter stub: module stub    
    """

    report_file = stub + '/'+  stub + '.html'
    param['module_report'] = open(param['working_dir']+'report/'+report_file, 'w')
    param['module_report'].write('<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 '+
                                'Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1'+
                                '-strict.dtd"><head><title></title></head><body>\n')
    param['module_report'].write('<center><h1>' + title + 'Overview</h1>')
    HELPER.write_html_table(param,
                                      table,
                                      out=param['module_report'],
                                      cell_width=80,
                                      fcol_width=150,
                                      deg=315)
    param['module_report'].write('<a href="' + stub +'_stats.txt">' + title +
                                ' statistics as tab delimited txt file</a>')
    #create an eSet:
    create_eset(out_file,
                param['pheno_file'],
                param,
                stub)
    HELPER.report_finish(param['module_report'])

    #add the fastqc html to the report
    param['report'].write('<a href="'+report_file+'">Full report</a><br>')

