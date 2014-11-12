#  Copyright (c) 2014, Boston University. All rights reserved.
#  
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are met: 
#  
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer. 
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution. 
#  
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
#  ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
#  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#  
#  The views and conclusions contained in the software and documentation are those
#  of the authors and should not be interpreted as representing official policies, 
#  either expressed or implied, of Boston University.
#  
#  Authors:
#    Daniel Gusenleitner [1,2], Vinay Kartha [1,2], Francesca Mulas [2], 
#    Yuxiang Tan [1,2], Liye Zhang [2], Stefano Monti [1,2]
#
#  [1] Bioinformatics Program, Boston University
#  [2] Center for Computational Biomedicine, Boston University  
#  import os, sys,  time, shlex, subprocess, module_helper, helper

def initialize_single_cpu(param):
    #split module list and add module load to the list
    helper.checkParameter(param,key='modules',dType=str)
    param['modules']=[('module load '+mod.strip()) for mod in param['modules'].split(",")]

        
def run_single_job(index, param, job_id, py_file, cores):
    #fetch modules that need to be loaded and add the python command for a single sample
    command_list = param['modules'][:]
    command_list.append('python2.7 '+ param['scripts_dir']+ py_file +' -i ' + str(index) + ' -n $NSLOTS')

    #call qsub script
    for cmd in command_list:
        args = shlex.split(cmd)
        output,error = subprocess.Popen(args,stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()

