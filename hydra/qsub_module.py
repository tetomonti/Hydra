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

""" qsub module that handles all interaction with the qsub system, which includes
writing a qsub file, submitting to the qsub system and waiting for the jobs to
finish.
"""

import os
import time
import subprocess
import hydra.module_helper
MODULE_HELPER = hydra.module_helper
from hydra.logs import writeLog

def initialize_qsub(param):
    """ Init function that initializes all qsub parameters

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    MODULE_HELPER.check_parameter(param, key='qsub_email', dtype=str)
    MODULE_HELPER.check_parameter(param, key='qsub_send_email', dtype=bool)
    MODULE_HELPER.check_parameter(param, key='qsub_memory', dtype=str)
    MODULE_HELPER.check_parameter(param, key='qsub_PROJECT', dtype=str)
    MODULE_HELPER.check_parameter(param, key='qsub_MACHINE', dtype=str)
    MODULE_HELPER.check_parameter(param, key='qsub_RUNTIME_LIMIT', dtype=str)
    MODULE_HELPER.check_parameter(param, key='qsub_wait_time', dtype=int)
    MODULE_HELPER.check_parameter(param, key='qsub_num_processors', dtype=str)


def wait_for_qsub(param, job_id):
    """This function checks the qsub system and determines whether the job has finished

    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    :Parameter job_id: random integer representing the job id for the current batch
    """

    user, _ = subprocess.Popen(args=['whoami'],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE).communicate()
    user = user.rstrip()

    #wait until there are no jobs with the current job name in the queue anymore
    qjobs = param['num_samples']
    writeLog('Waiting for single '+
             param['current_flag']+
             ' jobs to finish.... \n',
             param)
    while qjobs > 0:
        qjobs = 0
        time.sleep(param['qsub_wait_time'])
        # check the number of jobs
        call = ['qstat', '-u', '$USER']
        output, _ = subprocess.Popen(call,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE).communicate()

        #add grabbing the job IDs
        if len(output) == 0:
            qjobs = 0
        else:
            for current_line in output.split('\n'):
                if user in current_line and job_id in current_line:
                    qjobs += 1
    #change writing permissions in qsub directory so that group has access to
    #the qlog files. Otherwise other people cannot rerun the pipeline
    for filename in os.listdir(param['qsub_dir']):
        os.chmod(param['qsub_dir'] + filename, 0770)



def submit_jobs(index, param, py_file, job_id, cores, mem_free):
    """Function that submits a single job into the qsub system

    :Parameter index: index of the current file we are working on
    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    :Parameter py_file: wrapper that needs to be called to run the current step
    :Parameter cores: number of cores that should be used
    """

    cmd = (py_file +
           ' -i ' + str(index) +
           ' -n $NSLOTS' +
           ' -d ' +  param['working_dir'])

    #make a directory for all the qsub commands
    param['qsub_dir'] = param['working_dir']+'results/qsub/'
    if not os.path.exists(param['qsub_dir']):
        os.makedirs(param['qsub_dir'])

    qsub_dir = param['working_dir']+'results/qsub/'
    qsub_filename = qsub_dir+(param['stub'])[index]+'.qsub'

    outhandle = QsubClass(qsub_filename)
    outhandle.qsub_dir = qsub_dir
    outhandle.job_id = job_id
    outhandle.email = param['qsub_email']
    outhandle.send_email = param['qsub_send_email']
    if mem_free == 'standard':
        outhandle.memory = param['qsub_memory']
    else:
        outhandle.memory = mem_free
    outhandle.project = param['qsub_PROJECT']
    outhandle.machine = param['qsub_MACHINE']
    outhandle.runtime_limit = param['qsub_RUNTIME_LIMIT']
    outhandle.output_file([cmd])

    #call qsub script
    call = ['qsub']
    call.append('-pe')
    call.append('single_node')
    call.append(cores)
    call.append(qsub_filename)
    _, _ = subprocess.Popen(call,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE).communicate()

class QsubClass():
    """Class to output a qsub file

    """
    def __init__(self, name):
        """Standard initialization function.

        :Parameter name: filename of the qsub file
        """
        self.filename = name
        self.qsub_dir = './'
        self.job_id = -1
        self.email = ''
        self.machine = ''
        self.memory = ''
        self.runtime_limit = ''
        self.project = ''
        self.send_email = False


    def output_handle_gen(self, filepath=os.getcwd()):
        """Function that returns an output filehandle

        :Parameter filepath: complete path to the file that needs to be created
        """
        if '/' in self.filename:
            complete_path = self.filename
        else:
            complete_path = filepath + '/' + self.filename
        self.handle = open(complete_path, 'w')

    def output_file(self, command_line_list):
        """Function that writes to the qsub file

        :Parameter command_line_list: list of commands that should be executed
        """
        self.output_handle_gen()

        if self.machine == "scc":
            self.handle.write("source ~/.bashrc\n")
            #self.handle.write("#!bin/bash\n")
            self.handle.write("\n")

        else:
            self.handle.write("#!bin/bash\n")
            self.handle.write("#\n")
            self.handle.write("\n")

        if self.runtime_limit != '':
            self.handle.write("#$ -l h_rt="+self.runtime_limit+'\n')
            self.handle.write("\n")

        self.handle.write("#Specify which shell to use\n")
        self.handle.write("#$ -S /bin/bash\n")
        self.handle.write("\n")

        self.handle.write("#Run on the current working folder\n")
        self.handle.write("#$ -cwd\n")
        self.handle.write("\n")

        self.handle.write("#Give this job a name\n")
        self.handle.write("#$ -N "+self.job_id+'\n')
        self.handle.write("\n")

        self.handle.write("#Join standard output and error to a single file\n")
        self.handle.write("#$ -j y\n")
        self.handle.write("\n")

        self.handle.write("# Name the file where to redict standard output and error\n")
        if self.filename.count("/") >= 1:
            filename_info_list = self.filename.split("/")
            filename_info = filename_info_list[-1]
        else:
            filename_info = self.filename
        self.handle.write("#$ -o "+ self.qsub_dir + filename_info +".qlog\n")
        self.handle.write("\n")

        if self.project != "":
            self.handle.write("# Project this job belongs to \n")
            self.handle.write("#$ -P " + self.project+ " \n")
            self.handle.write("\n")

        if  self.email != "" and self.send_email:
            self.handle.write("# Send an email when the job begins and when it ends running\n")
            self.handle.write("#$ -m be\n")
            self.handle.write("\n")

            self.handle.write("# Whom to send the email to\n")
            self.handle.write("#$ -M "+self.email+ "\n")
            self.handle.write("\n")

        if self.memory != '':
            self.handle.write("# memory usage\n")
            self.handle.write("#$ -l mem_free="+self.memory+ "\n")
            self.handle.write("\n")

        self.handle.write("# Now let's Keep track of some'+ \
                          ' information just in case anything go wrong\n")
        self.handle.write("echo "+'"'+"========================================" + '"'+'\n')
        self.handle.write("echo "+'"'+"Starting on : $(date)"+'"'+ "\n")
        self.handle.write("echo "+'"'+"Running on node : $(hostname)"+'"'+"\n")
        self.handle.write("echo "+'"'+"Current directory : $(pwd)"+'"'+"\n")
        self.handle.write("echo "+'"'+"Current job ID : $JOB_ID"+'"'+"\n")
        self.handle.write("echo "+'"'+"Current job name : $JOB_NAME"+'"'+"\n")
        self.handle.write("echo "+'"'+"Task index number : $TASK_ID"+'"'+"\n")
        self.handle.write("echo "+'"'+"========================================" + '"'+'\n')
        self.handle.write("\n")

        for command_line in command_line_list:
            self.handle.write(command_line)
            self.handle.write('\n')

        self.handle.write("\n")
        self.handle.write("echo "+'"'+"========================================" + '"'+'\n')
        self.handle.write("echo "+'"'+"Finished on : $(date)"+'"'+ "\n")
        self.handle.write("echo "+'"'+"========================================" + '"'+'\n')
        self.handle.close()
