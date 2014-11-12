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
#  

import module_helper, os

def init(param):
   module_helper.checkParameter(param,key='match_pairs_exec',dType=str, checkFile=True)
   module_helper.checkParameter(param,key='match_pairs_python_version',dType=str)


def run_match_pairs(param,infile,infile2,outfile,outfile2):
    #run the script that matches pairs

    call = [param['match_pairs_python_version'],param['match_pairs_exec'],param[infile],param[infile2],outfile,outfile2]
    output,error = subprocess.Popen(call ,stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()

    param['file_handle'].write(error)
    param['file_handle'].write(output)

    with open(outfile+'.txt','w') as f:
        f.write(output)

    # Error handling
    if (not os.path.exists(outfile)) or os.stat(outfile).st_size<1000000:
        sys.exit(0)


if __name__ == "__main__":
    import subprocess, sys, os
    param=module_helper.initialize_module()

    #run match pairs
    outfile = param['module_dir']+param['stub'][param['file_index']]+'.clipped.matched.fastq.gz'
    outfile2 = param['module_dir']+param['stub'][param['file_index']]+'.clipped.matched.2.fastq.gz'

    run_match_pairs(param,'working_file','working_file2',outfile,outfile2)
    module_helper.wrapup_module(param,[outfile,outfile2])

