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

"""This module runs other modules such as tophat or cutadapt in serial instead of
parallelizing them.
"""
import shlex
import subprocess

def run_single_job(index, py_file):
    """Runs jobs sequentially on the current node

    :Parameter index: index of the current file
    :Parameter job_id: not used
    """

    #right now I'm using only one core in single cpu mode
    call = [py_file]
    call.append('-i')
    call.append(str(index))
    call.append('-n')
    call.append('1')
    
    args = shlex.split(' '.join(call))
    _, _ = subprocess.Popen(args,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE).communicate()

