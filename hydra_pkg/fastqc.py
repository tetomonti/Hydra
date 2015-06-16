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
from hydra_pkg import helper as HELPER
from hydra_pkg import module_helper as MODULE_HELPER



def remove_failed(param, input_files):
    #remove the failed samples from the run log
    failed = get_samples_to_copy(param, input_files)
    failed = [sample.split('/')[0] for sample in failed]
    
    if len(failed) > 0:
        HELPER.writeLog('The following samples did not pass the QC filter: '+
                        ' '.join(failed)+'\n',
                        param)
        for idx in range(len(param['stub'])):
            if param['stub'][idx] in failed:
                param['run_log'][-1][idx] = False

def get_faulty_samples(param, fqc_dir):
    #extract the names of samples that did not pass the sequence quality check
    overview = []
    for stub in param['fastqc_stub']:
        #extract checkmarks
        csv_file = open(fqc_dir+stub+'/summary.txt')
        csv_reader = csv.reader(csv_file, delimiter='\t')
        checks = [row[0] for row in csv_reader]
        csv_file.close()
        if checks[1] not in ['PASS']:
            overview.append(stub)
    return overview


def get_samples_to_copy(param, input_files):
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
    return(get_faulty_samples(param, fqc_dir))


def copy_files(param, input_files):
    """Copies all relevant fastqc files from the results directory into the report directory

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    #if there is no fastqc directory in the report make one
    fqc_dir = param['working_dir']+'results/fastqc/'
    param['fastqc_dir'] = param['working_dir']+'report/fastqc/'
    #Create a report fastqc subdirectory for each sample
    for stub in param['stub']:
        if not os.path.exists(param['fastqc_dir']+'/'+stub):
            os.makedirs(param['fastqc_dir']+'/'+stub)

    #depending on the flag either copy all finished samples or only the ones with issues
    if param['include_full_fastqc_report']:
        copy = param['fastqc_stub']
    else:
        copy = param['fastqc_overview']

    for fastqc_file in copy:
        call = ['cp', '-R']
        call.append(fqc_dir + fastqc_file+'/')
        call.append(param['fastqc_dir']+fastqc_file.split('/')[0])
        output, error = subprocess.Popen(call,
                                         stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE).communicate()
  

def create_overview_table(param, out, sub_mode=True):
    """Function that creates an html overview table continaing the most important fastqc statistics

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """

    if sub_mode:
        samples = param['fastqc_stub']
        icon_base = '../'
        link_base = ''
        include_report = param['include_full_fastqc_report']
    else:
        samples = param['fastqc_overview']
        icon_base = ''
        link_base = 'fastqc/'
        include_report = True

    #create a table
    table = []
    #put in headers
    table.append([stub.split('/')[0] for stub in samples])

    #link to summary files
    if include_report:
        temp = ['Summary files']
        for stub in samples:
            temp.append('<a href="'+link_base+stub+
                        '/fastqc_data.txt">raw</a>')
        table.append(temp)

        #link to overview files
        temp = ['Full report']
        for stub in samples:
            temp.append('<a href="'+link_base+stub+
                        '/fastqc_report.html"><img src="'+icon_base+'Icons/fastqc_icon.png"></a>')
        table.append(temp)

    #extract check marks
    table = table+extract_tables(param, icon_base, link_base, samples, include_report)
    #write the table as html
    HELPER.write_html_table(param, table, out, cell_width=30)


def extract_tables(param, icon_base, link_base, samples, include_report):
    """Extracts all relevant information for the overview table and writes it

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """

    fqc_dir = param['working_dir']+'results/fastqc/'
    #get the rownames
    csv_file = open(fqc_dir+samples[0]+'/summary.txt')
    csv_reader = csv.reader(csv_file, delimiter='\t')
    #get the rownames
    checks = [[row[1]] for row in csv_reader]
    csv_file.close()

    #links to the icons
    pass_icon = '<img src="'+icon_base+'Icons/tick.png">'
    fail_icon = '<img src="'+icon_base+'Icons/error.png">'
    warn_icon = '<img src="'+icon_base+'Icons/warning.png">'

    #get the values for each sample (icons for pass, faile or warning)
    for idx in range(len(samples)):
        csv_file = open(fqc_dir+samples[idx]+'/summary.txt')
        overview_file = link_base + samples[idx]+'/fastqc_report.html#M'
        csv_reader = csv.reader(csv_file, delimiter='\t')

        i = 0
        for row in csv_reader:
            if include_report:
                cell = '<a href="'+overview_file+str(i)+'">'
            else:
                cell = ''
            if row[0] == 'PASS':
                cell = cell+pass_icon
            if row[0] == 'FAIL':
                cell = cell+fail_icon
            if row[0] == 'WARN':
                cell = cell+warn_icon
            if include_report:
                cell = cell+'</a>'
            checks[i].append(cell)
            i += 1
        csv_file.close()
    return checks


def check_and_condense_list(file_list):
    confirmed = []
    for fil in file_list:
        if os.path.exists(fil):
            confirmed.append(fil)
    return confirmed

def read_raw_fastqc(param, output_files):
    """Reads raw fastqc outputs and processes it

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    :Parameter output_files: filenames to write the summaries in
    """
    summary_files = []
    for idx in range(len(param['fastqc_stub'])):
        summary_files.append(param['fastqc_dir']+
                             param['fastqc_stub'][idx]+
                             '/summary.txt')

    #check if the summary files actually exist
    summary_files = check_and_condense_list(summary_files)

    if len(summary_files) > 0:
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

        #depending on the flag either copy all finished samples or only the ones with issues
        if param['include_full_fastqc_report']:
            copy = param['fastqc_stub']
        else:
            copy = param['fastqc_overview']

        for idx in range(len(copy)):
            data_files.append(param['fastqc_dir']+
                              copy[idx]+
                              '/fastqc_data.txt')

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
        filehandle.write(' \t'+'\t'.join(copy)+'\n')
        for nam in key_list:
            filehandle.write(nam+'\t'+'\t'.join([str(vv) \
                             for vv in param['fast_qc_summary'][nam]])+'\n')
        filehandle.close()


def get_bar_width(fig_width, param):
    #calculates the bar plot width
    bar_width = 0.8 
    if len(param['fastqc_stub']) > 20:
        bar_width = fig_width / float(len(param['fastqc_stub'])) * 10
    return bar_width


def plot_number_of_reads(param, output_files, out):
    """Creates a plot with the number of reads per sample and plots it

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    :Parameter output_files: filenames to save th plot in
    """

    fqc_dir = param['working_dir']+'results/fastqc/'
    summary_files = []
    for idx in range(len(param['fastqc_stub'])):
        summary_files.append(fqc_dir+
                             param['fastqc_stub'][idx]+
                             '/fastqc_data.txt')
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
    fig_width = min (MODULE_HELPER.get_max_image_width(), 3+len(param['fastqc_stub'])*0.4)
    fig.set_size_inches(fig_width, 8)
    index = np.arange(len(num_total_reads))
    bar_width = get_bar_width(fig_width, param)
    opacity = 0.4
    _ = plt.bar(index,
                num_total_reads,
                bar_width,
                alpha=opacity,
                color='b')
    plt.xlabel('Samples')
    plt.ylabel('Number of reads')
    plt.title('Total number of reads')
    ticks = param['fastqc_stub']

    if fig_width != MODULE_HELPER.get_max_image_width():
        plt.xticks(index + bar_width / 2, ticks, rotation='vertical')
    plt.tight_layout()

    #put it into the report
    filename = 'report/fastqc/'+output_files+'_total_reads.png'
    pylab.savefig(param['working_dir']+filename)
    out.write('<img src="'+
              output_files+
              '_total_reads.png" '+
              'alt="total number of reads"><br><br>\n')

def plot_gc_content(param, input_files, out):
    """Creates a plot with the gc content per sample and plots it

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    :Parameter output_files: filenames to save th plot in
    """
    #extract number of reads, these are needed not only here but also for the bamqc
    summary_files = []
    fqc_dir = param['working_dir']+'results/fastqc/'
    for idx in range(len(param['fastqc_stub'])):
        summary_files.append(fqc_dir+
                             param['fastqc_stub'][idx]+
                             '/fastqc_data.txt')

    gc_content = []
    for sum_file in summary_files:
        gc_content.append(open(sum_file).readlines()[9].split('\t')[1].strip().rstrip())
    gc_content = [int(num) for num in gc_content]

    #create plot
    fig, axis = plt.subplots()
    fig_width = min (MODULE_HELPER.get_max_image_width(), 3+len(param['fastqc_stub'])*0.4)
    fig.set_size_inches(fig_width, 8)
    index = np.arange(len(gc_content))
    bar_width = get_bar_width(fig_width, param)
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
    if fig_width != MODULE_HELPER.get_max_image_width():
        plt.xticks(index + bar_width / 2, ticks, rotation='vertical')
    plt.tight_layout()
    axis.set_ylim(0, 100)

    #put it into the report
    filename = 'report/fastqc/'+input_files+'_gc_content.png'
    pylab.savefig(param['working_dir']+filename)
    out.write('<img src="'+
              input_files+
              '_gc_content.png" alt="GC content"><br><br>\n')


def write_sub_report(param, input_files, header):
    #separate report html for the fastqc results
    report_file = 'fastqc/'+input_files+'_fastqc.html'
    param['fastqc_report'] = open(param['working_dir']+'report/'+report_file, 'w')
    param['fastqc_report'].write('<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 '+
                                 'Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1'+
                                 '-strict.dtd"><head><title></title></head><body>\n')
    param['fastqc_report'].write('<center><h1>' + header + '</h1>')
    param['fastqc_report'].write('<a href="fastqc/'+
                                 input_files+
                                 'overview.txt">Table as tab delimited file'+
                                 '</a></center><br><br><br><br><br><br>')

    create_overview_table(param, param['fastqc_report'])
    plot_number_of_reads(param, input_files, param['fastqc_report'])
    plot_gc_content(param, input_files, param['fastqc_report'])

    HELPER.report_finish(param['fastqc_report'])
    return report_file


def report(param, input_files='fastq_files', header='FastQC results'):
    """Writes out html report with all relevant fastqc stats

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    :Parameter input_files: flag indicating which files to use (param key)
    :Parameter header: plain text that contains the header for the report
    """

    #assemble the full fastqc report
    param['fastqc_overview'] = get_samples_to_copy(param, input_files)
    copy_files(param, input_files)
    param['num_total_reads'] = [0]*param['num_samples']
    if len(param['fastqc_stub']) > 0:
        param['report'].write('<center><br><h2>'+header+'</h2><br><br><br>')
        read_raw_fastqc(param, input_files)

        #overview table only includes samples with issues      
        if len(param['fastqc_overview']) == 0:
            param['report'].write('All samples passed per base sequence quality check<br>')
        else:
            create_overview_table(param, param['report'], sub_mode=False)

        #write subreport
        report_file = write_sub_report(param, input_files, header)

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
    MODULE_HELPER.check_parameter(param, key='include_full_fastqc_report', dtype=bool)
    MODULE_HELPER.check_parameter(param, key='remove_failed', dtype=bool)


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


