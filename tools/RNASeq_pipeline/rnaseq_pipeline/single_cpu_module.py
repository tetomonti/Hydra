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
    cmd = py_file +' -i ' + str(index) + ' -n 1'
    args = shlex.split(cmd)
    _, _ = subprocess.Popen(args, stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE).communicate()

