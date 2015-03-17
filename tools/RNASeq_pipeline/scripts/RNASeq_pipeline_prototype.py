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
    from rnaseq_pipeline import pipeline, helper
    import subprocess
    
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
 
    
    
    
