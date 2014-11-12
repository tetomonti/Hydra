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

# -*- coding: cp936 -*-
"""
Usage: python RNASeq_pipeline.py -p parameter_file.txt
-h help
-p parameter file				*[No default value]
"""


if __name__ == "__main__":
## Import modules
    import matplotlib
    matplotlib.use('Agg')
    import sys, getopt
    import pipeline, helper, subprocess
    
    list_args = sys.argv[1:]
    
    if '-h' in list_args:
        print __doc__
        sys.exit(0)
    elif len(list_args) <2 or '-p' not in list_args:
        print __doc__
        sys.exit(0)
    else:
        param, parameter_file, updates = helper.update_parameters(list_args)
       
    
    #if we resume we specify the old parameter file
    helper.initialize_standard(param)
    helper.write_updated_file(updates, param, parameter_file)

    helper.initialize_qsub(param)
    helper.readFastqFilenames(param)
    helper.initialize_logfiles(param)
    
    helper.writeLog('Initializing all module parameters ... \n',param)
    pipeline.initialize_all(param) 

    helper.writeLog('Initializing successful!\n',param)    
    helper.writeLog('####################################################\n',param)    
    
## run pipeline
    helper.writeLog('Running all modules: \n\n',param)
    pipeline.run_all(param)
   
## reporting
    helper.report_start(param) 
    pipeline.report_all(param)
    helper.report_finish(param)
 
    
    
    