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

"""Fastqc module
This module contains functions for running fastqc as
well as multiple reporting and plotting functions to add the QC into the
report html
"""
import os
import subprocess
import sys
import csv
import numpy as np
import matplotlib.pyplot as plt
import pylab
import hydra.helper
import hydra.module_helper
MODULE_HELPER = hydra.module_helper



def copy_files(param, input_files):
    """Copies all relevant fastqc files from the results directory into the report directory

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    #if there is no fastqc directory in the report make one
    param['fastqc_dir'] = param['working_dir']+'report/fastqc/'
    #Create a report fastqc subdirectory for each sample
    for stub in param['stub']:
        if not os.path.exists(param['fastqc_dir']+'/'+stub):
            os.makedirs(param['fastqc_dir']+'/'+stub)

    #reconstruct file names
    param['fastqc_stub'] = []

    for stb, filename in zip(param['stub'], param[input_files]):
        param['fastqc_stub'].append(stb+'/'+filename.split('/')[-1])

    if param['paired']:
        for stb, filename in zip(param['stub'], param[input_files+'2']):
            param['fastqc_stub'].append(stb+'/'+filename.split('/')[-1])

    #fix the suffix
    for idx in range(len(param['fastqc_stub'])):
        param['fastqc_stub'][idx] = param['fastqc_stub'][idx].replace('.fastq.gz',
                                                                      '_fastqc')

    #use only the files that exist
    fqc_dir = param['working_dir']+'results/fastqc/'
    param['fastqc_stub'] = [fn for fn in param['fastqc_stub'] if os.path.exists(fqc_dir+fn)]

    #copy the unpacked directories
    for fastqc_file in param['fastqc_stub']:
        call = ['cp', '-R']
        call.append(fqc_dir + fastqc_file+'/')
        call.append(param['fastqc_dir']+fastqc_file.split('/')[0])
        output, error = subprocess.Popen(call,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE).communicate()
    #cut off the suffix
    param['fastqc_stub'] = [fn.replace('_fastqc', '') for fn in param['fastqc_stub']]



def create_overview_table(param, out):
    """Function that creates an html overview table continaing the most important fastqc statistics

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    #create a table
    table = []
    #put in headers
    table.append([stub for stub in param['fastqc_stub']])
    #link to summary files

    temp = ['Summary files']
    for stub in param['fastqc_stub']:
        temp.append('<a href="'+stub+
                    '_fastqc/fastqc_data.txt">raw</a>')
    table.append(temp)

    #link to overview files
    temp = ['Full report']
    for stub in param['fastqc_stub']:
        temp.append('<a href="'+stub+
                    '_fastqc/fastqc_report.html"><img src="../Icons/fastqc_icon.png"></a>')
    table.append(temp)
    #extract check marks
    table = table+extract_tables(param)
    #write the table as html
    MODULE_HELPER.write_html_table(param, table, out)

def extract_tables(param):
    """Extracts all relevant information for the overview table and writes it

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    #get the rownames
    csv_file = open(param['fastqc_dir']+param['fastqc_stub'][0]+'_fastqc/summary.txt')
    csv_reader = csv.reader(csv_file, delimiter='\t')
    #get the rownames
    checks = [[row[1]] for row in csv_reader]
    csv_file.close()

    #links to the icons
    pass_icon = '<img src="../Icons/tick.png">'
    fail_icon = '<img src="../Icons/error.png">'
    warn_icon = '<img src="../Icons/warning.png">'

    #get the values for each sample (icons for pass, faile or warning)
    for idx in range(len(param['fastqc_stub'])):
        csv_file = open(param['fastqc_dir']+param['fastqc_stub'][idx]+'_fastqc/summary.txt')
        overview_file = 'fastqc/'+param['fastqc_stub'][idx]+'_fastqc/fastqc_report.html#M'
        csv_reader = csv.reader(csv_file, delimiter='\t')

        i = 0
        for row in csv_reader:
            cell = '<a href="'+overview_file+str(i)+'">'
            if row[0] == 'PASS':
                cell = cell+pass_icon
            if row[0] == 'FAIL':
                cell = cell+fail_icon
            if row[0] == 'WARN':
                cell = cell+warn_icon
            cell = cell+'</a>'
            checks[i].append(cell)
            i += 1
        csv_file.close()
    return checks


def read_raw_fastqc(param, output_files):
    """Reads raw fastqc outputs and processes it

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    :Parameter output_files: filenames to write teh summaries in
    """
    summary_files = []
    for idx in range(len(param['fastqc_stub'])):
        summary_files.append(param['fastqc_dir']+
                             param['fastqc_stub'][idx]+
                             '_fastqc'+
                             '/summary.txt')

    fastqc = dict()
    #add entries into fastqc dictionary
    filehandle = open(summary_files[0])
    fastqc = dict([(name.split('\t')[1].strip(), []) for name in filehandle.readlines()])
    filehandle.close()

    #fill fastqc dictionary with information from the summary files
    for sum_file in summary_files:
        filehandle = open(sum_file)
        for name in filehandle.readlines():
            fastqc[name.split('\t')[1].rstrip()].append(name.split('\t')[0].strip())
        filehandle.close()

    key_list = fastqc.keys()

    #fill fastqc dictionary with information from the data file
    data_files = []
    for idx in range(len(param['fastqc_stub'])):
        data_files.append(param['fastqc_dir']+
                          param['fastqc_stub'][idx]+
                          '_fastqc/fastqc_data.txt')

    labels = ['Encoding',
              'Total Sequences',
              'Filtered Sequences',
              'Sequence length',
              '%GC']

    fastqc.update(dict([(l, []) for l in labels]))
    key_list.extend(labels)

    for d_file in data_files:
        filehandle = open(d_file)
        for name in filehandle.readlines():
            if name.split('\t')[0].strip() in fastqc.keys():
                fastqc[name.split('\t')[0].strip()].append(name.split('\t')[1].rstrip())
        filehandle.close()

    param['fast_qc_summary'] = fastqc

    # write overview file
    filehandle = open(param['fastqc_dir']+output_files+'overview.txt', 'w')
    filehandle.write(' \t'+'\t'.join(param['fastqc_stub'])+'\n')
    for nam in key_list:
        filehandle.write(nam+'\t'+'\t'.join([str(vv) \
                         for vv in param['fast_qc_summary'][nam]])+'\n')
    filehandle.close()

def plot_number_of_reads(param, output_files, out):
    """Creates a plot with the number of reads per sample and plots it

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    :Parameter output_files: filenames to save th plot in
    """
    summary_files = []
    for idx in range(len(param['fastqc_stub'])):
        summary_files.append(param['fastqc_dir']+
                             param['fastqc_stub'][idx]+
                             '_fastqc/fastqc_data.txt')
    #print summary_files
    num_total_reads = []
    for sum_file in summary_files:
        num_total_reads.append(open(sum_file).readlines()[6].split('\t')[1].strip().rstrip())
    num_total_reads = [int(num) for num in num_total_reads]

    if not param['paired']:
        param['num_total_reads'] = num_total_reads
    else:
        param['num_total_reads'] = [0]*param['num_samples']

        for idx in range(param['num_samples']):
            param['num_total_reads'][idx] = (num_total_reads[idx]+
                                             num_total_reads[idx+param['num_samples']])

    #create plot
    fig, _ = plt.subplots()
    fig.set_size_inches(3+len(param['fastqc_stub'])*0.4, 8)
    index = np.arange(len(num_total_reads))
    bar_width = 0.8
    opacity = 0.4
    _ = plt.bar(index,
                num_total_reads,
                bar_width,
                alpha=opacity,
                color='b')
    plt.xlabel('Samples')
    plt.ylabel('Total number of reads')
    plt.title('Total number of reads across samples')
    ticks = param['fastqc_stub']
    plt.xticks(index + bar_width, ticks, rotation='vertical')
    plt.tight_layout()

    #put it into the report
    filename = 'report/fastqc/'+output_files+'_total_reads.png'
    pylab.savefig(param['working_dir']+filename)
    out.write('<img src="'+
              output_files+
              'total_reads.png" '+
              'alt="total number of reads"><br><br>\n')

def plot_gc_content(param, input_files, out):
    """Creates a plot with the gc content per sample and plots it

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    :Parameter output_files: filenames to save th plot in
    """
    #extract number of reads, these are needed not only here but also for the bamqc
    summary_files = []
    for idx in range(len(param['fastqc_stub'])):
        summary_files.append(param['fastqc_dir']+
                             param['fastqc_stub'][idx]+
                             '_fastqc/fastqc_data.txt')

    gc_content = []
    for sum_file in summary_files:
        gc_content.append(open(sum_file).readlines()[9].split('\t')[1].strip().rstrip())
    gc_content = [int(num) for num in gc_content]

    #create plot
    fig, axis = plt.subplots()
    fig.set_size_inches(3+len(param['fastqc_stub'])*0.4, 8)
    index = np.arange(len(gc_content))
    bar_width = 0.8
    opacity = 0.4
    _ = plt.bar(index,
                gc_content,
                bar_width,
                alpha=opacity,
                color='b')
    plt.xlabel('Samples')
    plt.ylabel('%GC content')
    plt.title('GC content across samples')
    ticks = param['fastqc_stub']
    plt.xticks(index + bar_width, ticks, rotation='vertical')
    plt.tight_layout()
    axis.set_ylim(0, 100)

    #put it into the report
    filename = 'report/fastqc/'+input_files+'_gc_content.png'
    pylab.savefig(param['working_dir']+filename)
    out.write('<img src="'+
              input_files+
              'gc_content.png" alt="GC content"><br><br>\n')


def report(param, input_files='fastq_files', header='FastQC results'):
    """Writes out html report with all relevant fastqc stats

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    :Parameter input_files: flag indicating which files to use (param key)
    :Parameter header: plain text that contains the header for the report
    """

    #assemble the full fastqc report
    copy_files(param, input_files)
    param['num_total_reads'] = [0]*param['num_samples']
    if len(param['fastqc_stub']) > 0:
        param['report'].write('<center><br><h2>'+header+'</h2>')
        read_raw_fastqc(param, input_files)


        #separate report html for the fastqc results
        report_file = 'fastqc/'+input_files+'_fastqc.html'
        param['fastqc_report'] = open(param['working_dir']+'report/'+report_file, 'w')
        param['fastqc_report'].write('<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 '+
                                     'Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1'+
                                     '-strict.dtd"><head><title></title></head><body>\n')
        param['fastqc_report'].write('<center><h1>' + header + '</h1></center>')

        create_overview_table(param, param['fastqc_report'])
        param['fastqc_report'].write('<a href="fastqc/'+
                                     input_files+
                                     'overview.txt">Table as tab delimited file'+
                                     '</a><br><br><br>')
        plot_number_of_reads(param, input_files, param['fastqc_report'])
        plot_gc_content(param, input_files, param['fastqc_report'])


        hydra.helper.report_finish(param['fastqc_report'])

        #add the fastqc html to the report
        param['report'].write('<a href="'+report_file+'">Full report</a><br>')

    else:
        param['report'].write('There were no results to show.')


def init(param):
    """Initialization function, that checks if the bamqc_script that is run
    on every single samples is available

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    MODULE_HELPER.check_parameter(param, key='fastqc_exec', dtype=str)

def run_fastqc(filename, param):
    """Wrapper that build the command line call

    :Parameter filename: param key pointing to the inputfile list
    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    out_dir = param['module_dir']+'/'+param['outstub']
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    #build command line call
    call = [param['fastqc_exec']]
    call.append(param[filename])
    call.append('-o')
    call.append(out_dir)
    call.append('--extract')
    param['file_handle'].write('CALL: '+' '.join(call)+'\n')

    output, error = subprocess.Popen(call,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE).communicate()
    param['file_handle'].write(error)
    param['file_handle'].write(output)
    if 'Analysis complete for' not in output:
        sys.exit()


def main():
    """Main function that is run on each samples, which in turn calls the
    actual fastqc tool to extract the QC statistics
    """
    param = MODULE_HELPER.initialize_module()

    #run fastqc
    run_fastqc('working_file', param)

    if not param['paired']:
        MODULE_HELPER.wrapup_module(param)

    #calling it on the second fastq file if it is paired
    else:
        run_fastqc('working_file2', param)
        MODULE_HELPER.wrapup_module(param)


