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

"""featureCount module
This module contains functions for initializing all featureCount specific variables and
a wrapper that that runs tophat using those parameters on a single sample.
In addition it also contains functions to extract and write statistics and
a wrapper that calls an R script
"""

from hydra_pkg import module_helper as MODULE_HELPER
from hydra_pkg import helper as HELPER
import os
import re
import subprocess



def init(param):
    """Initialization function that checks the all relevant tophat parameters

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    MODULE_HELPER.check_parameter(param, key='featureCount_exec', dtype=str)
    MODULE_HELPER.check_parameter(param, key='featureCount_s', dtype=str)
    MODULE_HELPER.check_parameter(param, key='featureCount_t', dtype=str)
    MODULE_HELPER.check_parameter(param, key='featureCount_id', dtype=str)
    MODULE_HELPER.check_parameter(param, key='featureCount_gft', dtype=str, checkfile=True)
    MODULE_HELPER.check_parameter(param, key='featureCount_by_meta', dtype=bool)
    MODULE_HELPER.check_parameter(param, key='Rscript_exec', dtype=str)

def process_stat_files(param):
    """Copies all relevant files into the report directory and also extracts
    the total number of reads from the bamqc output

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """

    #get the files that are actually in the output directory
    call = ['cp', '-R']
    call.append(param['working_dir']+'results/featureCount/')
    call.append(param['working_dir']+'report/')
    _, _ = subprocess.Popen(call,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE).communicate()

    featurecount_file = (param['working_dir']+
                         'results/featureCount/featureCount_stats.txt')
    #extract table
    table = []
    filehandle = open(featurecount_file)
    #header
    table.append(filehandle.readlines()[0].rstrip().split('\t'))
    table[0] = table[0][1:]
    filehandle.close()

    #total number of aligned reads
    tot_reads = param['bam_qc']['unique_aligned_reads']

    filehandle = open(featurecount_file)
    for line in filehandle.readlines()[1:]:
        cur_line = line.rstrip().split('\t')
        cur_line[0] = re.sub(r'_',' ',cur_line[0])
        if cur_line[0] != 'Unassigned MultiMapping':
            perc = ([cur_line[0]]+
                    MODULE_HELPER.get_percentage(cur_line[1:],
                                                 tot_reads,
                                                 len(cur_line)-1))
        table.append(perc)
    filehandle.close()
    return table


  
def report(param):
    """Function that writes all HTSeq related statistics into the html report

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    featurecount_dir = param['working_dir']+'report/featureCount/'
    if not os.path.exists(featurecount_dir):
        os.makedirs(featurecount_dir)
        
    #report only if there were actually results
    out_file = param['working_dir']+'deliverables/featureCount_raw_counts.txt'
    if os.path.exists(out_file):
        param['report'].write('<center><br><br><h2>FeatureCount statistics</h2>')
        table = process_stat_files(param)
        MODULE_HELPER.create_sub_report(param, out_file, table, 'featureCount', 'FeatureCount')                                
        MODULE_HELPER.plot_count_overview(param, 'featureCount', table)


def finalize(param, input_files='count_files'):
    """This function is run after featureCount is run on each sample. It collects all results
    and puts them into a file

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    :Parameter input_files: flag that indicates the input files
    """

    import csv
    HELPER.writeLog('Collecting featureCount raw counts ... \n', param)

    #check which of these files are actually available
    working_files = [iFile for iFile in param[input_files] if iFile != '']

    if len(working_files) > 0:
        #get feature ID using the first column in the first file in the list of working files
        csv_file = open(working_files[0])
        csv_reader = csv.reader(csv_file, delimiter='\t')

        #For featureCount output, we want to skip the first two lines as they
         #include the featureCount call and the headers which we don't want
        next(csv_reader, None)
        next(csv_reader, None)

        #Now start by taking the list of identifier,
        #which is the first column in the file
        counts = [row[0] for row in csv_reader]
        csv_file.close()

        #get all the expression values
        header = 'ID'
        for idx in range(param['num_samples']):
            if param[input_files] != '':
                header = header+'\t'+param['stub'][idx]
                csv_file = open(param[input_files][idx])
                csv_reader = csv.reader(csv_file, delimiter='\t')

                #Here too we want to skip the first two lines, before getting the counts
                next(csv_reader, None)
                next(csv_reader, None)
                #Now start getting the counts (row[6]) and add in the ID (counts[i]) before it
                idx = 0
                for row in csv_reader:
                    counts[idx] = counts[idx]+'\t'+row[6]
                    idx += 1
                csv_file.close()

        #output the file
        out_file = param['working_dir']+'deliverables/featureCount_raw_counts.txt'
        out_handle = open(out_file, 'w')
        out_handle.write(header+'\n')

        for i in range(len(counts)):
            out_handle.write(counts[i]+'\n')
        out_handle.close()

        #output_phenotype_file
        HELPER.writeLog('Writing phenotype data ... \n', param)
        MODULE_HELPER.output_sample_info(param)

        #write summary stats
        #featureCount does this on its own so we can just fetch each summary file
        #check which of these files are actually available
        working_files = [iFile+'.summary' for iFile in param[input_files] if iFile != '']

        if len(working_files) > 0:
            #get Status column from summary file using the first column in
            #the first file in the list of working files
            csv_file = open(working_files[0])
            csv_reader = csv.reader(csv_file, delimiter='\t')
            #Here, we want to skip the first line, as it simply points to the
            #alignment file used when running featureCount
            next(csv_reader, None)

            #Now start by taking the list of identifier,
            #which is the first column in the file
            entry = [row[0] for row in csv_reader]
            csv_file.close()

            #get all the summary stats for each sample
            header = 'Status'

            for idx in range(param['num_samples']):
                if param[input_files] != '':
                    header = header+'\t'+param['stub'][idx]
                    #Fetch the corresponding sample's summary file
                    csv_file = open(param[input_files][idx]+'.summary')
                    csv_reader = csv.reader(csv_file, delimiter='\t')
                    #Again, we want to skip the first line, as it simply points
                    #to the alignment file used when running featureCount
                    next(csv_reader, None)

                    #Now start getting the stats (row[1]) and add in the Status
                    # (counts[i]) before it
                    i = 0
                    for row in csv_reader:
                        entry[i] = entry[i]+'\t'+row[1]
                        i += 1
                    csv_file.close()
            #output the file
            out_handle = open(param['working_dir']+
                              'results/featureCount/featureCount_stats.txt',
                              'w')
            out_handle.write(header+'\n')

            for i in range(len(entry)):
                out_handle.write(entry[i]+'\n')
            out_handle.close()
        else:
            print 'featureCount was not run successfully on any of the files..\n'



def main():
    """Main function that is run on each samples, which in turn calls runs
    featureCount on a sample.
    """
    import sys
    param = MODULE_HELPER.initialize_module()


    outfile = param['module_dir']+param['outstub']
    call = [param['featureCount_exec']]


    if param['paired']:
        call = call + ['-P', '-p']
    if param['featureCount_by_meta'] == False:
        call.append('-f')

    call = call + ['-t', param['featureCount_t']]
    call = call + ['-s', param['featureCount_s']]
    call = call + ['-g', param['featureCount_id']]
    call = call + ['-a', param['featureCount_gft']]
    call = call + ['-o', outfile]
    call.append(param['working_file'])

    param['file_handle'].write('CALL: '+' '.join(call)+'\n')
    output, error = subprocess.Popen(call,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE).communicate()

    param['file_handle'].write(error)
    param['file_handle'].write(output)

    if not os.path.exists(outfile):
        param['file_handle'].write('featureCount run failed \n')
        sys.exit(0)

    #wrap up and return the current workingfile
    MODULE_HELPER.wrapup_module(param, [outfile])



