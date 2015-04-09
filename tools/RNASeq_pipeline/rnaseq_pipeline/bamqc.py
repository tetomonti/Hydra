"""Bamqc module
This module contains functions for running the custom bamqc script as
well as multiple reporting and plotting functions to add the QC into the
report html
"""
import os
import subprocess
import rnaseq_pipeline.module_helper
MODULE_HELPER = rnaseq_pipeline.module_helper
import numpy as np
import matplotlib.pyplot as plt
import pylab

def copy_files(param):
    """Copies all relevant bamqc files from the results directory into the report directory

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    #if there is no bamqc directory in the report make one
    param['bamqc_dir'] = param['working_dir'] + 'report/bamqc/'
    if not os.path.exists(param['bamqc_dir']):
        os.makedirs(param['bamqc_dir'])

    #get the files that are actually in the output directory
    call = 'ls ' + param['working_dir'] + 'results/bamqc/'
    output, _ = subprocess.Popen(call.split(),
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE).communicate()
    param['bamqc_stub'] = [line for line in output.split('\n') if line != '']

    #use only the files that are in the stub
    param['bamqc_stub'] = [bqc_stub \
                          for bqc_stub in param['bamqc_stub'] \
                          if bqc_stub in param['stub']]

    #copy the unpacked directories
    for stub in param['bamqc_stub']:
        call = 'cp -R ' + \
               param['working_dir'] + \
               'results/bamqc/' + \
               stub + ' ' + \
               param['bamqc_dir']
        output, _ = subprocess.Popen(call.split(),
                                         stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE).communicate()

def create_overview_table(param):
    """Function that creates an html overview table continaing the most important bamQC statistics

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    #create a table
    table = []
    table.append([stub for stub in param['bamqc_stub']])
    num_s = len(param['bamqc_stub'])
    #link to summary files
    table.append(['Summary files'] + \
                 ['<a href="bamqc/' + \
                 stub + \
                 '/output.txt">raw</a>' for stub in param['bamqc_stub']])
    #link to overview files
    table.append(['Full report'] + \
             ['<a href="bamqc/' + \
             stub + \
             '/sample_stats.html"><img src="Icons/fastqc_icon.png"></a></td>' \
             for stub in param['bamqc_stub']])
    #header
    table.append(['Percentages based on total number of reads:'])
    #percent aligned
    table.append(['Percent aligned'] + \
                 MODULE_HELPER.getPercentage(\
                 number1=param['bam_qc']['single_count_alignments'],
                 number2=param['num_total_reads'],
                 ntotal=num_s))
    #percent uniquely aligned
    table.append(['Percent uniquely aligned'] + \
                 MODULE_HELPER.getPercentage(\
                 number1=param['bam_qc']['unique_aligned_reads'],
                 number2=param['num_total_reads'],
                 ntotal=num_s))
    #header
    table.append(['Percentages based on total number of alignments:'])
    #percent single reads
    table.append(['Percent single end reads'] + \
                 MODULE_HELPER.getPercentage(\
                 number1=param['bam_qc']['is_singleton'],
                 number2=param['bam_qc']['total_aligned_reads'],
                 ntotal=num_s))
    #percent paired reads
    table.append(['Percent paired end reads'] + \
                 MODULE_HELPER.getPercentage(\
                 number1=param['bam_qc']['is_paired'],
                 number2=param['bam_qc']['total_aligned_reads'],
                 ntotal=num_s))
    #percent proper paired reads
    table.append(['Percent proper paired reads'] + \
                 MODULE_HELPER.getPercentage(\
                 number1=param['bam_qc']['is_proper_pair'],
                 number2=param['bam_qc']['total_aligned_reads'],
                 ntotal=num_s))
    #percent spliced reads
    table.append(['Percent spliced reads'] + \
                 MODULE_HELPER.getPercentage(\
                 number1=param['bam_qc']['spliced_reads'],
                 number2=param['bam_qc']['total_aligned_reads'],
                 ntotal=num_s))
    #percent insert reads
    table.append(['Percent of reads with inserts'] + \
                 MODULE_HELPER.getPercentage(\
                 number1=param['bam_qc']['reads_with_inserts'],
                 number2=param['bam_qc']['total_aligned_reads'],
                 ntotal=num_s))
    #percent deletion reads
    table.append(['Percent of reads with deletions'] + \
                 MODULE_HELPER.getPercentage(\
                 number1=param['bam_qc']['reads_with_deletions'],
                 number2=param['bam_qc']['total_aligned_reads'],
                 ntotal=num_s))
    MODULE_HELPER.writeHTMLtable(param, \
                                 table, \
                                 out=param['report'], \
                                 cell_width=65)

def read_raw_bamqc(param):
    """Reads the raw output from the bamqc run of a single bamqc run

    :Parameter param: - dictionary that contains all general RNASeq pipeline parameters
    """

    summary_files = [param['bamqc_dir'] + \
                     param['bamqc_stub'][idx] + \
                     '/output.txt' for idx in range(len(param['bamqc_stub']))]
    bamqc = dict()
    #add entries into bamqc dictionary
    filehandle = open(summary_files[0])
    name_list = []
    for name in filehandle.readlines():
        bamqc[name.split('\t')[0].strip()] = []
        name_list.append(name.split('\t')[0].strip())
    filehandle.close()

    #fill bamqc dictionary
    for sum_file in summary_files:
        filehandle = open(sum_file)
        for name in filehandle.readlines():
            bamqc[name.split('\t')[0].strip()].append(\
            float(name.split('\t')[1].rstrip()))
        filehandle.close()

    #get the actual counts
    param['bam_qc'] = bamqc
    param['bam_qc']['single_count_alignments'] = \
           [sum([v[i] for k, v in bamqc.items() \
           if 'num_multiread' in k]) \
           for i in range(len(bamqc['num_multiread 1']))]
    param['bam_qc']['num_total_reads'] = param['num_total_reads']

    key_list = ['num_total_reads', 'single_count_alignments']
    key_list.extend(name_list)

    filehandle = open(param['bamqc_dir'] + 'overview.txt', 'w')
    filehandle.write(' \t' + '\t'.join(param['bamqc_stub']) + '\n')
    for nam in key_list:
        filehandle.write(nam + '\t' + \
                 '\t'.join([str(vv) for vv in param['bam_qc'][nam]])+'\n')
    filehandle.close()

def plot_alignment_overview(param):
    """Creates a plot that contains the statistic on the number of aligned reads

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    #values:
    unique = param['bam_qc']['unique_aligned_reads'][:]
    aligned = param['bam_qc']['single_count_alignments'][:]
    total = param['bam_qc']['num_total_reads'][:]
    total = [(total[i]-aligned[i]) for i in range(len(param['bamqc_stub']))]
    aligned = [(aligned[i]-unique[i]) for i in range(len(param['bamqc_stub']))]

    #create plot
    fig, _ = plt.subplots()
    fig.set_size_inches(5 + len(param['bamqc_stub']) * 0.4, 8)
    index = np.arange(len(param['bamqc_stub']))

    bar_width = 0.8
    opacity = 0.4
    rects1 = plt.bar(index,
                     total,
                     bar_width,
                     bottom=param['bam_qc']['single_count_alignments'],
                     alpha=opacity,
                     color='b')
    rects2 = plt.bar(index,
                     aligned,
                     bar_width,
                     bottom=unique,
                     alpha=opacity,
                     color='r')
    rects3 = plt.bar(index,
                     unique,
                     bar_width,
                     alpha=opacity,
                     color='g')

    plt.xlabel('Samples')
    plt.ylabel('Total (aligned) reads')
    plt.title('Number of reads across samples')
    ticks = param['bamqc_stub']
    plt.xticks(index + bar_width, ticks, rotation='vertical')
    plt.legend((rects1[0],
                rects2[0],
                rects3[0]),
                ('Total reads', 'Aligned reads', 'Uniquely aligned'), \
                loc='lower left')
    plt.tight_layout()

    #put it into the report
    filename = 'report/bamqc/aligned_reads.png'
    pylab.savefig(param['working_dir'] + filename)
    param['report'].write('<img src="bamqc/aligned_reads.png"' + \
                          ' alt="number of aligned reads"><br><br>\n')

def plot_spliced_reads(param):
    """Creates a plot that contains the statistic on the number of spliced reads

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """

    percent = [0.0]*len(param['bamqc_stub'])
    for idx in range(len(param['bamqc_stub'])):
        percent[idx] = round(float(param['bam_qc']['spliced_reads'][idx]) /  \
                       float(param['bam_qc']['single_count_alignments'][idx])\
                       * 100, 3)

    #create plot
    fig, _ = plt.subplots()
    fig.set_size_inches(5 + len(param['bamqc_stub']) * 0.4, 8)
    index = np.arange(len(param['bamqc_stub']))

    bar_width = 0.8
    opacity = 0.4
    _ = plt.bar(index, percent, bar_width, alpha=opacity, color='b')
    plt.xlabel('Samples')
    plt.ylabel('Percentage of spliced reads of all aligned reads')
    plt.title('Percentage of spliced reads across samples')
    ticks = param['bamqc_stub']
    plt.xticks(index + bar_width, ticks, rotation='vertical')
    plt.tight_layout()

    #put it into the report
    filename = 'report/bamqc/spliced_reads.png'
    pylab.savefig(param['working_dir'] + filename)
    param['report'].write('<img src="bamqc/spliced_reads.png" ' + \
                          'alt="number of spliced reads"><br><br>\n')


def plot_insert_reads(param):
    """Creates a plot that contains the statistic on the number of reads with inserts

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    percent = [0.0]*len(param['bamqc_stub'])
    for idx in range(len(param['bamqc_stub'])):
        percent[idx] = round(\
                       float(param['bam_qc']['reads_with_inserts'][idx])/ \
                       float(param['bam_qc']['single_count_alignments'][idx]) \
                       * 100, 3)

    #create plot
    fig, _ = plt.subplots()
    fig.set_size_inches(5 + len(param['bamqc_stub']) * 0.4, 8)
    index = np.arange(len(param['bamqc_stub']))

    bar_width = 0.8
    opacity = 0.4
    _ = plt.bar(index, percent, bar_width, alpha=opacity, color='b')
    plt.xlabel('Samples')
    plt.ylabel('Percent of reads with inserst of all aligned reads')
    plt.title('Percent of reads with inserts across samples')
    ticks = param['bamqc_stub']
    plt.xticks(index + bar_width, ticks, rotation='vertical')
    plt.tight_layout()

    #put it into the report
    filename = 'report/bamqc/insert_reads.png'
    pylab.savefig(param['working_dir']+filename)
    param['report'].write('<img src="bamqc/insert_reads.png" '+ \
                          'alt="number of inserted reads"><br><br>\n')

def plot_delete_reads(param):
    """Creates a plot that contains the statistic on the number of reads that contain deletions

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    percent = [0.0] * len(param['bamqc_stub'])
    for idx in range(len(param['bamqc_stub'])):
        percent[idx] = round( \
                       float(param['bam_qc']['reads_with_deletions'][idx]) / \
                       float(param['bam_qc']['single_count_alignments'][idx]) \
                       * 100, 3)

    #create plot
    fig, _ = plt.subplots()
    fig.set_size_inches(5 + len(param['bamqc_stub']) * 0.4, 8)
    index = np.arange(len(param['bamqc_stub']))

    bar_width = 0.8
    opacity = 0.4
    _ = plt.bar(index, percent, bar_width, alpha=opacity, color='b')
    plt.xlabel('Samples')
    plt.ylabel('Percent of reads with deletions of all aligned reads')
    plt.title('Percent of reads with deletion across samples')
    ticks = param['bamqc_stub']
    plt.xticks(index + bar_width, ticks, rotation='vertical')
    plt.tight_layout()

    #put it into the report
    filename = 'report/bamqc/delete_reads.png'
    pylab.savefig(param['working_dir']+filename)
    param['report'].write('<img src="bamqc/delete_reads.png"' + \
                          ' alt="number of deleted reads"><br><br>\n')

def plot_paired_singleton(param):
    """Creates a plot that contains the statistics on the number of paired, \
    proper paired and singleton reads

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    #values:
    single = param['bam_qc']['is_singleton'][:]
    paired = param['bam_qc']['is_paired'][:]
    proper = param['bam_qc']['is_proper_pair'][:]
    paired_bottom = [(single[i] + proper[i]) \
                    for i in range(len(param['bamqc_stub']))]
    paired = [paired[i]-proper[i] for i in range(len(param['bamqc_stub']))]

    #create plot
    fig, _ = plt.subplots()
    fig.set_size_inches(5 + len(param['bamqc_stub']) * 0.4, 8)
    index = np.arange(len(param['bamqc_stub']))

    bar_width = 0.8
    opacity = 0.4
    rects1 = plt.bar(index, paired, bar_width, \
                     bottom=paired_bottom, alpha=opacity, color='b')
    rects2 = plt.bar(index, proper, bar_width, \
                     bottom=single, alpha=opacity, color='r')
    rects3 = plt.bar(index, single, bar_width, alpha=opacity, color='g')

    plt.xlabel('Samples')
    plt.ylabel('Number of single / paired reads')
    plt.title('Number of single/paired/proper paired reads across samples')
    ticks = param['bamqc_stub']
    plt.xticks(index + bar_width, ticks, rotation='vertical')
    plt.legend((rects1[0], \
                rects2[0], \
                rects3[0]), \
               ('Paired end reads', \
                'Proper paired reads', 'Single end reads'), loc='lower left')
    plt.tight_layout()

    #put it into the report
    filename = 'report/bamqc/paired_reads.png'
    pylab.savefig(param['working_dir']+filename)
    param['report'].write('<img src="bamqc/paired_reads.png"' + \
                          ' alt="number of paired reads"><br><br>\n')

def plot_mismatches(param):
    """Creates a plot that split the reads by number of mismatches

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    #values:
    mm0 = [param['bam_qc']['num_unique_mismatches 0'][i]+param['bam_qc']\
          ['num_multiple_mismatches 0'][i]\
          for i in range(len(param['bamqc_stub']))]
    mm1 = [param['bam_qc']['num_unique_mismatches 1'][i]+param['bam_qc']\
          ['num_multiple_mismatches 1'][i]\
          for i in range(len(param['bamqc_stub']))]
    mm2 = [param['bam_qc']['num_unique_mismatches 2'][i]+param['bam_qc']\
          ['num_multiple_mismatches 2'][i]\
          for i in range(len(param['bamqc_stub']))]
    mm3 = [param['bam_qc']['num_unique_mismatches 3'][i]+param['bam_qc']\
          ['num_multiple_mismatches 3'][i]\
          for i in range(len(param['bamqc_stub']))]
    mm4 = [param['bam_qc']['num_unique_mismatches 4'][i]+param['bam_qc']\
          ['num_multiple_mismatches 4'][i]\
          for i in range(len(param['bamqc_stub']))]

    #make it cummulative
    bot_mm2 = [mm1[i]+mm0[i] for i in range(len(param['bamqc_stub']))]
    bot_mm3 = [bot_mm2[i]+mm2[i] for i in range(len(param['bamqc_stub']))]
    bot_mm4 = [bot_mm3[i]+mm3[i] for i in range(len(param['bamqc_stub']))]

    #create plot
    fig, _ = plt.subplots()
    fig.set_size_inches(5 + len(param['bamqc_stub']) * 0.4, 8)
    index = np.arange(len(param['bamqc_stub']))

    bar_width = 0.8
    opacity = 0.4
    rects1 = plt.bar(index, mm0, bar_width,
                     alpha=opacity, color='b')
    rects2 = plt.bar(index, mm1, bar_width,
                     bottom=mm0, alpha=opacity, color='r')
    rects3 = plt.bar(index, mm2, bar_width,
                     bottom=bot_mm2, alpha=opacity, color='g')
    rects4 = plt.bar(index, mm3, bar_width,
                     bottom=bot_mm3, alpha=opacity, color='#555555')
    rects5 = plt.bar(index, mm4, bar_width, bottom=bot_mm4,
                  alpha=opacity, color='#ffff00')

    plt.xlabel('Samples')
    plt.ylabel('Mismatches')
    plt.title('Number of mismatches across samples')
    ticks = param['bamqc_stub']
    plt.xticks(index + bar_width, ticks, rotation='vertical')
    plt.legend((rects1[0],
                rects2[0],
                rects3[0],
                rects4[0],
                rects5[0]),
                ('Perfect match',
                 '1 mismatch',
                 '2 mismatches',
                 '3 mismatches',
                 '4+ mismatches'), loc='lower left')
    plt.tight_layout()

    #put it into the report
    filename = 'report/bamqc/mismatches.png'
    pylab.savefig(param['working_dir']+filename)
    param['report'].write('<img src="bamqc/mismatches.png" ' + \
                          'alt="number of mismatches"><br><br>\n')

def report(param):
    """This function creates a full html report for the bamQC and also \
    copies all relevant files

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    param['report'].write('<center><br><br><br><br><h2>Bam QC results</h2>')
    copy_files(param)
    if len(param['bamqc_stub']) > 0:
        read_raw_bamqc(param)
        create_overview_table(param)
        param['report'].write('<a href="bamqc/overview.txt">' + \
                     'QC results as tab delimited file</a><br><br><br>')
        plot_alignment_overview(param)
        plot_mismatches(param)
        plot_paired_singleton(param)
        plot_spliced_reads(param)
        plot_insert_reads(param)
        plot_delete_reads(param)
    else:
        param['report'].write('There were no results to show.')


def init(param):
    """Initialization function, that checks if the bamqc_script that is run \
    on every single samples is available

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    MODULE_HELPER.checkParameter(param, key='bamqc_script', dType=str)


def main():
    """Main function that is run on each samples, which in turn calls the \
    actual bamqc running script to extract the QC statistics
    """
    import sys
    param = MODULE_HELPER.initialize_module()

    #run create output directory
    outdir = param['module_dir']+param['stub'][param['file_index']]+'/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)


    call = param['bamqc_script'] + ' -i ' + \
                    param['working_file'] + ' -o ' + outdir

    param['file_handle'].write('CALL: '+call+'\n')
    output, error = subprocess.Popen(call.split(),
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE).communicate()
    param['file_handle'].write(error)
    param['file_handle'].write(output)

    #if outputfile doesn't exits?
    if not os.path.exists(outdir+'stats.json'):
        param['file_handle'].write('QC did not finish correctly..')
        sys.exit()

    MODULE_HELPER.wrapup_module(param)

