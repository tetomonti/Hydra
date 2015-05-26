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

import os
import sys

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

    if optional and key not in param:
        param[key] = ''
    else:
        #check if key is present
        if key not in param:
            print 'Parameter '+key+' is missing in the parameter file.'
            sys.exit(0)

        #cast to correct data type
        if dtype == bool:
            param[key] = param[key] in ['True', 'TRUE', 'true', 'T', '1']
        else:
            param[key] = dtype(param[key])

        #check if the value for the key is allowed
        if len(allowed) > 0 and param[key] not in allowed:
            print 'Parameter '+key+' can only be one of the following: ' + allowed
            sys.exit(0)

        #if file or directory check if it exists
        if checkfile and not os.path.exists(param[key]):
            print 'The file in '+key+' = ', param[key], ' does not exist'
            sys.exit(0)


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

    #Import modules
    import getopt
    import json

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


def rotate_word(word, deg=270):
    """Creates the html code to display a word rotated (by default 270 degrees).
    Of note is that this works in all standard browsers

    :Parameter word: the word that should be printed
    :Parameter deg: the degree of the rotation
    """
    return('<div style="float: center;position: relative;-moz-transform: '+
           'rotate(270deg);  /* FF3.5+ */-o-transform: rotate(270deg);'+
           '  /* Opera 10.5 */ -webkit-transform: rotate('+str(deg)+
           'deg);  /* Saf3.1+, Chrome */ filter:  progid:DXImageTransform.'+
           'Microsoft.BasicImage(rotation=3);  /* IE6,IE7 */ -ms-filter: '+
           'progid:DXImageTransform.Microsoft.BasicImage(rotation=3);' +
           ' /* IE8 */">' +word+'</div>')

def write_html_table(param, table, out, fcol_width=200, cell_width=50, initial_breaks=8, deg=315):
    """HTML table writing function that takes a table and writes it nicely out
    as an html table

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    :Parameter table: table to output
    :Parameter out: output file handle
    :Parameter fcol_width: width of the first column
    :Parameter cell_width: width of each cell, except first column
    :Parameter initial_breaks: spaces before the column
    :Parameter deg: rotation of the text in the header line
    """

    out.write('<table id="one-column-emphasis" class="fixed"><col width="'+
              str(fcol_width)+'px"/>\n')
    out.write(''.join(['<br>']*initial_breaks)+'\n')
    out.write(''.join(['<col width="'+str(cell_width)+'px"/>\n']*len(table[0])))
    #write headertable
    out.write('<thead><tr><th></th>'+
              ''.join(['<th>'+
              rotate_word(stub.replace('-', '_'), deg)+
                          '</th>\n' for stub in table[0]])+
              '</tr></thead>\n')

    #write the pass and fail for each module
    for idx in range(1, len(table)):
        out.write('<tr>')
        if len(table[idx]) == 1:
            out.write('<th colspan="'+str(1+len(table[0]))+
                      '">'+table[idx][0]+'</th>\n')
        else:
            for i in range(len(table[idx])):
                out.write('<td>'+str(table[idx][i])+'</td>\n')
        out.write('</tr>\n')

    #close table
    out.write('</table><br>')

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
