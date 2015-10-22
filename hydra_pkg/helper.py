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

"""pipeline helper functions. Mixed bag of reading, writing, submission
and checking functions all related to general pipeline functionality
"""
import json
import os
import random
import re
import shutil
import subprocess
import sys
from hydra_pkg.logs import writeLog
from hydra_pkg import qsub_module
from hydra_pkg import single_cpu_module
from hydra_pkg import RNASEQ_PIPELINE_DIR



###############################################################################
#########       Init Functions  ###############################################
###############################################################################

def check_parameter(param, key, dtype, allowed=[], checkfile=False, optional=False):
    """Parameter checking wrapper which calls that module helper function

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

def initialize_logfiles(param):
    """Function that creates the log files

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    #create log files and log directory if they do not exist already
    log_dir = param['working_dir']+'results/log/'
    #if the doesn't exist create it
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    #if log files don't exit create them
    for idx in range(param['num_samples']):
        log_file = log_dir+param['stub'][idx]+'.log'
        if not os.path.exists(log_file):
            open(log_file, 'a').close()
    param['log_handle'] = param['working_dir']+'results/main.log'

def initialize_qsub(param):
    """Function to initialize the parallelization module parameters

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    check_parameter(param, key='genome_annotation_gft', dtype=str, checkfile=True)    
    check_parameter(param, key='genome', dtype=str, checkfile=True)
    check_parameter(param, key='stranded', allowed=['no', 'yes', 'reversed'], dtype=str)    
    check_parameter(param, key='qsub_email', dtype=str)
    check_parameter(param, key='qsub_send_email', dtype=bool)
    check_parameter(param, key='qsub_memory', dtype=str)
    check_parameter(param, key='qsub_PROJECT', dtype=str)
    check_parameter(param, key='qsub_MACHINE', dtype=str)
    check_parameter(param, key='qsub_RUNTIME_LIMIT', dtype=str)
    check_parameter(param, key='qsub_wait_time', dtype=int)
    check_parameter(param, key='qsub_num_processors', dtype=str)

def clean_up(param):
    """Remove all old results

    :Parameter param: parameter object
    """
    #this is only relevant if the directories actually exist
    res_dir = os.path.exists(param['working_dir']+'results/')
    rep_dir = os.path.exists(param['working_dir']+'report/')
    del_dir = os.path.exists(param['working_dir']+'delivarables/')
    
    if res_dir | rep_dir | del_dir:
        #check before deleting any previous results
        if param['ask_before_deleting']:
            answer = raw_input('Are you sure you want to delete '+
                               'all existing results? (yes/no): ')
        else:
            answer = 'yes'

        if answer == 'yes':
            if os.path.exists(param['working_dir']+'results/'):
                shutil.rmtree(param['working_dir']+'results/')
            if os.path.exists(param['working_dir']+'report/'):
                shutil.rmtree(param['working_dir']+'report/')
            if os.path.exists(param['working_dir']+'deliverables/'):
                shutil.rmtree(param['working_dir']+'deliverables/')
        else:
            print ('Stopping... If you want to resume, '+
                   'please adjust the parameter file.')
            exit(0)


def clean_failed(param):
    """Remove log files of the samples that failed in the last step, this will
    force them to be rerun from scratch

    :Parameter param: parameter object
    """
    mainlog = param['working_dir']+'results/main.log'
    if os.path.exists(mainlog):
        filehandle = open(mainlog, 'r')
        error_samples = ''
        for line in filehandle:
            line = line.strip()
            if 'error in samples' in line:
                error_samples = line
        if error_samples is not '':
            error_samples = re.sub('error in samples ','',error_samples)
            error_samples = error_samples.split(';')
            #delete all log files
            for sample in error_samples:
                logfile = param['working_dir']+'results/log/'+sample+'.log'
                if os.path.exists(logfile):
                    os.remove(logfile)
            #and finally delete the mainlog
            os.remove(mainlog)


def initialize_standard(param):
    """checks default pipeline parameters and create working directories

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    check_parameter(param, key='raw_filenames', dtype=str, checkfile=True)
    check_parameter(param, key='paired', dtype=bool)
    check_parameter(param, key='clean_run', dtype=bool)
    check_parameter(param, key='ask_before_deleting', dtype=bool)
    check_parameter(param, key='remove_intermediate', dtype=bool)
    check_parameter(param, key='run_failed_from_scratch', dtype=bool)
    check_parameter(param, key='verbose', dtype=bool)
    check_parameter(param, key='run_single_cpu', dtype=bool)
    check_parameter(param, key='raw_file_header', dtype=bool)
    check_parameter(param, key='aligner', dtype=str)
    check_parameter(param, key='QC_and_trim_only', dtype=bool)
    check_parameter(param, key='run_cufflinks', dtype=bool)
    check_parameter(param, key='run_htseq', dtype=bool)
    check_parameter(param, key='run_featureCount', dtype=bool)
    check_parameter(param, key='zipped_fastq', dtype=bool)

    qsub_module.initialize_qsub(param)

    #checking working directory and going there
    check_parameter(param, key='working_dir', dtype=str, checkfile=True)
    if  param['working_dir'][len(param['working_dir'])-1] != '/':
        param['working_dir'] = param['working_dir']+'/'

    #if directory exists and the pipeline
    #should be run from scratch delete the directory
    if param['clean_run']:
        clean_up(param)

    if param['run_failed_from_scratch']:
        clean_failed(param)

    #if results or report directory do not exist create them
    if not os.path.exists(param['working_dir']+'results/'):
        os.makedirs(param['working_dir']+'results/')
    if not os.path.exists(param['working_dir']+'report/'):
        os.makedirs(param['working_dir']+'report/')
    if not os.path.exists(param['working_dir']+'deliverables/'):
        os.makedirs(param['working_dir']+'deliverables/')


def parse_parameters(par_file):
    """This function constructs a dictionary with pairs: parameter name - value

    :Parameter par_file: parameter file location
    """

    dict_param = dict([(l.split("=")[0].strip(),
                        l.split("#")[0].split("=")[1].strip()) \
                        for l in par_file.readlines() \
                        if "=" in l and l.strip()[0] != "#"])
    return dict_param

def read_fastq_filenames(param):
    """Get filenames and stubs, while checking if the files actually exist in
    paired mode there will be 2 fastq files per sample, otherwise only 1

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    param['stub'] = []
    param['raw_files'] = []
    if param['paired']:
        param['raw_files2'] = []

    #assign raw file locations
    filehandle = open(param['raw_filenames'], 'r')
    #skip first line if there is a header specified
    header = param['raw_file_header']
    for line in filehandle:
        if header:
            header = False
        else:
            line = line.strip()
            line = line.split('\t')
            param['raw_files'].append(line[0])
            # check if filename actually exists
            if not os.path.exists(line[0]):
                print 'The file '+line[0]+' does not exist'
                sys.exit(0)
            if not param['paired']:
                param['stub'].append(line[1])
            else:
                param['raw_files2'].append(line[1])
                #check if file exists
                if not os.path.exists(line[0]):
                    print 'The file '+line[1]+' does not exist'
                    sys.exit(0)
                param['stub'].append(line[2])
    param['num_samples'] = len(param['stub'])

    #start a log file that keeps track on which files are successfully completed
    param['run_log'] = [[True]*param['num_samples']]
    param['run_log_headers'] = ['raw']

def update_parameters(args):
    """This function updates the parameter dictionaries with
    additional values provided in the command line

    :Parameter args: arguments provided in the command line
    """
    pnames = args[0::2]
    tuple_args = zip(pnames, args[1::2])
    file_param = filter(lambda x: x[0] == '-p', tuple_args)
    # read parameter file
    fname = file_param[0][1]
    param = parse_parameters(open(fname))
    # read additional parameters
    other_param = filter(lambda x: x[0] != '-p', tuple_args)
    new_param = []
    if len(other_param) > 0:
        # check additional parameter names
        for t_arg in other_param:
            if '--' in t_arg[0] and t_arg[0].split('--')[1] in param.keys():
                new_param.append((t_arg[0].split('--')[1], t_arg[1]))
            else:
                print "Parameter %s is specified incorrectly, " \
                  "please see parameter file for all possible parameters." %t_arg[0]
                sys.exit(0)
    # update parameters
    new_dict = dict([t for t in new_param if t[0] in param.keys()])
    param.update(new_dict)
    return param, fname, new_dict

def write_updated_file(updated, param, parameter_file):
    """Writes out the updated parameter file

    :Parameter param: updated parameter dictionary
    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    :Parameter parameter_file: parameter filename

    """
    new_parameter_file = (param['working_dir']+
                          'results/'+
                          (parameter_file.split('/')[-1])[:-4]+
                          '_used.txt')
    with open(parameter_file) as infile:
        with open(new_parameter_file, "w") as outfile:
            for line in infile:
                k = line.split("=")[0].strip()
                if "=" in line and k in updated.keys():
                    outfile.write('%s= %s\n' %(line.split("=")[0], updated[k]))
                else:
                    outfile.write(line)
    param['parameter_file'] = new_parameter_file

###############################################################################
#####       Running Functions  ###############################################
###############################################################################

def add_input_output_files(param, input_files, output_files):
    """Store information on which file type to work and what output file to use

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    :Parameter input_files: input filenames
    :Parameter output_files: output filenames
    """

    #store information on which file type to work and what output file to use
    param['input_files'] = input_files
    param['output_files'] = output_files

    #check if the specified input files actually exist otherwise throw an error
    if not param.has_key(param['input_files']):
        writeLog('The input file type you specified does not exist', param)
        sys.exit(0)

    #check if the key for the specified output files actually exist otherwise create them
    if not param.has_key(param['output_files']) and param['output_files'] != '':
        param[param['output_files']] = ['']*param['num_samples']

def check_queue_success(param):
    """Check how many jobs ran successful and prints that into the main log

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    :return parameter: failed number of samples
    """
    s_noend = []
    for idx in range(param['num_samples']):
        logfile = open(param['working_dir']+'results/log/'+param['stub'][idx]+'.log')
        lines_end = [line for line in logfile.readlines() \
                      if 'ENDING %s |' %(param['current_flag']) in line.rstrip()]
        #check if the job was successfully finished
        if len(lines_end) == 1:
            #check if there is suppose to be a output file and if that is
            #the case store the location of the output file
            if  param['output_files'] != '':
                retval = lines_end[0].split('|')[1].strip()
                if len(retval.split(';')) == 1:
                    param[param['output_files']][idx] = retval
                else:
                    param[param['output_files']][idx] = retval.split(';')[0]
                    param[param['output_files']+'2'][idx] = retval.split(';')[1]
            #add an entry into the run log that shows that the
            #current job has been finished successfully
            param['run_log'][-1][idx] = True
        else:
            if param['output_files'] != '':
                param[param['output_files']][idx] = ''
            s_noend.append(param['stub'][idx])

    #output number of unsucessfully run samples
    if len(s_noend) > 0:
        writeLog('error in samples %s' %(';'.join([s for s in s_noend])), param)
        writeLog('\n', param)
    else:
        writeLog(param['current_flag']+' successful!\n\n', param)
    return len(s_noend)

def dump_parameters(param):
    """dumps the parameter object into a JSON object

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    with open(param['working_dir']+'results/parameters.json', 'w') as filehandle:
        json.dump(param, filehandle)


def is_module_finished(param):
    """check if the current module was already run successfully before

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    if param['clean_run']:
        finished = False
    else:
        handle = open(param['log_handle'])
        success = [line for line in handle.readlines() \
                    if '%s successful!' %(param['current_flag']) in line.rstrip()]
        handle.close()
        finished = len(success) == 1
    return finished

def check_job_finished(param, index):
    """Function that checks if the current job is already finished

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    :Parameter index: index of the current file we are working on
    """   
    #check if it is not a clean run and if module already completed
    finished = False
    if not param['clean_run']:
        log_file = (param['working_dir']+'results/log/'+
                    param['stub'][index]+'.log')        
        file_handle = open(log_file)
        for line in file_handle.readlines():
            if 'ENDING %s |' %(param['current_flag']) in line.rstrip():
                finished = True
        file_handle.close()
    return finished


def submit_job(param, py_file, input_files, output_files='', cores='1-8', mem_free='standard'):
    """submit a job, which runs the same wrapper on every single file

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    :Parameter py_file: name of the python wrapper
    :Parameter input_files: input filenames
    :Parameter output_files: output filenames
    :Parameter cores: number of cores that should be used if available
    """

    #store information on which files to use as input and which ones as output
    add_input_output_files(param, input_files, output_files)

    #get the current flag for the log file
    # this is a terrible hard coded solution assuming that py_file has a value
    # 'run_current_dir' and we want to drop the 'run_' to get 'current_dir'
    param['current_dir'] = py_file[4:]
    param['current_flag'] = param['current_dir'] + '_' + input_files

    if is_module_finished(param):
        #add coment into main log file
        writeLog('Skipping '+
                 param['current_flag']+
                 ' since it was already successfully finished.\n',
                 param)
        #and fetch the current working files
        for idx in range(param['num_samples']):
            logfile = open(param['working_dir']+
                           'results/log/'+
                           param['stub'][idx]+
                           '.log')            
            #check if the current module actually produces output files that are
            #used subsequently. (i.e. tophat or the trimmer as opposed to fastqc or bamqc)
            if  param['output_files'] != '':
                #if there are output files fetch the values and add them to the
                #proper list of files (bam_files, fastq_files, ... )
                lines_end = [line for line in logfile.readlines() \
                             if 'ENDING %s | ' %(param['current_flag']) in line.rstrip()]
                retval = lines_end[0].split('|')[1].strip()
                #check if there are one or 2 working files returned.
                #e.g. cutadapt will return 2 values, while most modules return only 1
                if len(retval.split(';')) == 1:
                    param[param['output_files']][idx] = retval
                else:
                    param[param['output_files']][idx] = retval.split(';')[0]
                    param[param['output_files']+'2'][idx] = retval.split(';')[1]
                logfile.close()
        #add in the run log that all samples finished successfully
        param['run_log'].append([True]*param['num_samples'])
        param['run_log_headers'].append(param['current_flag'])

    else:
        writeLog('Starting '+param['current_flag']+'\n', param)
        #create current working directory
        param['module_dir'] = (param['working_dir']+
                               'results/'+
                               param['current_dir']+'/')
        if not os.path.exists(param['module_dir']):
            os.makedirs(param['module_dir'])
        #write all parameters to file so the single subnodes can use them
        dump_parameters(param)

        #generate a job id
        job_id = 'HyDrA_'+str(random.randint(1, 10000))

        #submit all jobs
        submitted = False
        for index in range(param['num_samples']):
            #first check if the job was already finished and skip if the pipeline is in resume mode
            #and also check if the step before was run sucessfully
            if not check_job_finished(param, index) and param['run_log'][-1][index]:
                if param['run_single_cpu']:
                    single_cpu_module.run_single_job(index, py_file)
                else:
                    qsub_module.submit_jobs(index,
                                            param,
                                            py_file,
                                            job_id,
                                            cores,
                                            mem_free)
                    submitted = True

        #wait for qsub scripts to finish if there were actually submitted jobs
        if submitted:
            qsub_module.wait_for_qsub(param, job_id)
        #start a new log column for the current batch of jobs
        param['run_log'].append([False]*param['num_samples'])
        param['run_log_headers'].append(param['current_flag'])
        #check if all jobs finished successful
        if check_queue_success(param) == param['num_samples']:
            writeLog('All samples failed, aborting... \n',param)            
            sys.exit(0)
        writeLog('++++++++++++++++++++++++++++\n\n', param)



###############################################################################
####       Reporting Functions  ###############################################
###############################################################################

def copy_and_link(param, param_key, text):
    """copies a file and creates and html link for it

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    call = ['cp', param[param_key], param['working_dir']+'report/']
    _, _ = subprocess.Popen(call,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE).communicate()
    #add file to html
    param['report'].write('<a href="'+
                          param[param_key].split('/')[-1]+
                          '">'+text+'</a><br>')



def report_finish(outhandle):
    """Writes the bottom of the html report

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    outhandle.write('<style>table.fixed { table-layout:fixed; } '+
                    'table.fixed td { overflow: hidden; }#one-column'+
                    '-emphasis{font-family:"Lucida Sans Unicode", '+
                    '"Lucida Grande", Sans-Serif;font-size:12px;width:'+
                    '480px;text-align:left;border-collapse:collapse;'+
                    'margin:20px;}#one-column-emphasis th{font-size:'+
                    '14px;font-weight:normal;color:#039;padding:12px'+
                    ' 15px;}#one-column-emphasis td{color:#669;border'+
                    '-top:1px solid #e8edff;padding:5px 5px;}.oce-'+
                    'first{background:#d0dafd;border-right:10px solid'+
                    ' transparent;border-left:10px solid transparent;}'+
                    '#one-column-emphasis tr:hover td{color:#339;'+
                    'background:#eff2ff;}</style></body>\n')
    outhandle.close()

def report_run_log(param):
    """Creates the run log html for the final report indicating which modules
    were run successfully on which samples.

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    param['report'].write('<center><br><h2>Running Log</h2>')

    #create a dedicated directory for the run-log
    log_dir = param['working_dir']+'report/run_log/'
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    #start the report
    report_file = 'run_log/run_log.html'
    param['runlog_report'] = open(param['working_dir']+'report/'+report_file, 'w')
    param['runlog_report'].write('<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 '+
                                 'Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1'+
                                 '-strict.dtd"><head><title></title></head><body>\n')
    param['runlog_report'].write('<center><h1>Run Log</h1></center>')

    #create a table
    table = []
    table.append([stub for stub in param['stub']])

    #links to the icons
    pass_icon = '<img src="../Icons/tick.png">'
    fail_icon = '<img src="../Icons/error.png">'

    #create table
    for idx in range(len(param['run_log_headers'])):
        line = [param['run_log_headers'][idx].title()]
        for i in range(param['num_samples']):
            status = fail_icon
            if param['run_log'][idx][i]:
                status = pass_icon
            line = line + [status]
        table.append(line)

    #write the table as html
    write_html_table(param, table, out=param['runlog_report'])
    report_finish(param['runlog_report'])

    #add the fastqc html to the report
    param['report'].write('<a href="'+report_file+'">Full report</a><br>')


def report_start(param):
    """Writes the header of the html report

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    writeLog('Starting to write the report ... \n', param)
    param['report'] = open(param['working_dir']+'report/index.html', 'w')
    param['report'].write('<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 '+
                          'Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1'+
                          '-strict.dtd"><head><title></title></head><body>\n')
    param['report'].write('<center><br><h1>Hydra Report</h1>')
    copy_and_link(param, 'parameter_file', 'Parameter file')
    copy_and_link(param, 'raw_filenames', 'Raw files')

    #copy the pass/fail icons into the directory
    call = ['cp', '-R']
    call.append(os.path.join(RNASEQ_PIPELINE_DIR, 'Icons'))
    call.append(param['working_dir']+'report/')
    _, _ = subprocess.Popen(call,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE).communicate()

    #report the run log in a table and show which samples passed/failed
    report_run_log(param)

def rotate_word(word, deg='315'):
    """Creates the html code to display a word rotated (by default 270 degrees).
    Of note is that this works in all standard browsers

    :Parameter word: the word that should be printed
    :Parameter deg: the degree of the rotation
    """
    return('<div style="float: center;position: relative;  bottom: 10px;'+
           '-ms-transform: rotate('+str(deg)+'deg); -webkit-transform:'+
           'rotate('+str(deg)+'deg);-moz-transform:  rotate('+str(deg)+'deg);-o-transform: '+
           'rotate('+str(deg)+'deg);writing-mode: rl-tb;/">' +str(word)+'</div>')


def write_html_table(param, table, out, fcol_width=250, cell_width=50, initial_breaks=5, deg=315):
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

