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

