import os, sys, shutil, json, time, shlex, subprocess, random 
import qsub_module, single_cpu_module, module_helper
  
################################################################################################################################################################################
################################################################################################################################################################################
####       Init Functions  ################################################################################################################################################
################################################################################################################################################################################
################################################################################################################################################################################
def checkParameter(param,key,dType,allowed=[],checkFile=False,optional=False):
    module_helper.checkParameter(param=param,key=key,dType=dType,allowed=allowed,checkFile=checkFile,optional=optional)
            
def initialize_logfiles(param):
    #create log files and log directory if they do not exist already
    log_dir=param['working_dir']+'results/log/'
    #if the doesn't exist create it
    if not os.path.exists(log_dir):
         os.makedirs(log_dir)
    #if log files don't exit create them
    for idx in range(param['num_samples']):
        log_file = log_dir+param['stub'][idx]+'.log'
        if not os.path.exists(log_file):
            open(log_file,'a').close()
    param['log_handle'] = param['working_dir']+'results/main.log'
        
def initialize_qsub(param):
    qsub_module.initialize_qsub(param)    

def initialize_standard(param):
    #check default pipeline parameters and create working directories
    checkParameter(param,key='raw_filenames',dType=str,checkFile=True)
    checkParameter(param,key='paired',dType=bool)
    checkParameter(param,key='clean_run',dType=bool)
    checkParameter(param,key='verbose',dType=bool)
    checkParameter(param,key='run_single_cpu',dType=bool)
    checkParameter(param,key='raw_file_header',dType=bool)
    checkParameter(param,key='aligner',dType=str)
    checkParameter(param,key='QC_and_trim_only',dType=bool)
    checkParameter(param,key='run_cufflinks',dType=bool)
    checkParameter(param,key='run_htseq',dType=bool)
    checkParameter(param,key='run_featureCount',dType=bool)
    
    #checking working directory and going there
    checkParameter(param,key='working_dir',dType=str,checkFile=True) 
    if  param['working_dir'][len(param['working_dir'])-1]!='/':
        param['working_dir']=param['working_dir']+'/'
  #  os.chdir(os.path.dirname(param['working_dir']))

    #directory where all the scripts are located
    checkParameter(param,key='scripts_dir',dType=str)
    if  param['scripts_dir'][len(param['scripts_dir'])-1]!='/':
        param['scripts_dir']=param['scripts_dir']+'/'

    #if directory exists and the pipeline should be run from scratch delete the directory
    if param['clean_run']:
        #check before deleting any previous results
        answer = input('Are you sure you want to delete all existing results? (yes/no): ')
        if answer != 'yes':
            print 'Stopping... If you want to resume, please adjust the parameter file.'
            exit(0)
        if os.path.exists(param['working_dir']+'results/'):
            shutil.rmtree(param['working_dir']+'results/')
        if os.path.exists(param['working_dir']+'report/'):
            shutil.rmtree(param['working_dir']+'report/')
        if os.path.exists(param['working_dir']+'deliverables/'):
            shutil.rmtree(param['working_dir']+'deliverables/')

    #if results or report directory do not exist create them
    if not os.path.exists(param['working_dir']+'results/'):
         os.makedirs(param['working_dir']+'results/')
    if not os.path.exists(param['working_dir']+'report/'):
         os.makedirs(param['working_dir']+'report/')
    if not os.path.exists(param['working_dir']+'deliverables/'):
         os.makedirs(param['working_dir']+'deliverables/')


#This function constructs a dictionary with pairs: parameter name - value
def parse_parameters(par_file):
    dict_param = dict([(l.split("=")[0].strip(), l.split("#")[0].split("=")[1].strip()) for l in par_file.readlines() if "=" in l and l.strip()[0]!="#"])
    return dict_param

def readFastqFilenames(param):
#Get filenames and stubs, while checking if the files actually exist in paired mode there will be 2 fastq files per sample, otherwise only 1
    param['stub']=[]
    param['raw_files']=[]
    if param['paired']: param['raw_files2']=[]
    
    #assign raw file locations
    f = open (param['raw_filenames'] , 'r')
    #skip first line if there is a header specified
    header=param['raw_file_header']
    for line in f:
        if header: 
            header=False
        else:
            line = line.strip()
            line = line.split('\t')
            param['raw_files'].append(line[0])
            # check if filename actually exists
            if not os.path.exists(line[0]):
                print 'The file '+line[0]+' does not exist';
                sys.exit(0)
            if not param['paired']:
                param['stub'].append(line[1])
            else:
                param['raw_files2'].append(line[1])
                #check if file exists
                if not os.path.exists(line[0]):
                    print 'The file '+line[1]+' does not exist';
                    sys.exit(0)
                param['stub'].append(line[2])
    param['num_samples']=len(param['stub'])
            
    #start a log file that keeps track on which files are successfully completed
    param['run_log']=[[True]*param['num_samples']]
    param['run_log_headers']=['raw']

def update_parameters(args):
    #This function updates the parameter dictionaries with additional values provided in the command line
    pnames = args[0::2]
    tuple_args = zip(pnames, args[1::2])
    file_param = filter(lambda x: x[0]== '-p', tuple_args)
    # read parameter file
    fname = file_param[0][1]
    param = parse_parameters(open(fname))
    # read additional parameters
    other_param = filter(lambda x: x[0]!= '-p', tuple_args)
    new_param = []
    if len(other_param) > 0:
        # check additional parameter names
        for t_arg in other_param:
            if '--' in t_arg[0] and t_arg[0].split('--')[1] in param.keys():
                new_param.append((t_arg[0].split('--')[1], t_arg[1]))
            else:
                print "Parameter %s is specified incorrectly, please see parameter file for all possible parameters." %t_arg[0]
                sys.exit(0) 
    # update parameters
    new_dict = dict([t for t in new_param if t[0] in param.keys()])
    param.update(new_dict)
    return param, fname, new_dict

def write_updated_file(updated, param, parameter_file):
    new_parameter_file=param['working_dir']+'results/'+(parameter_file.split('/')[-1])[:-4]+'_used.txt'
    with open(parameter_file) as f:
        with open(new_parameter_file, "w") as f1:
            for l in f:
                k = l.split("=")[0].strip()
                if "=" in l and k in updated.keys():
                    f1.write('%s= %s\n' %(l.split("=")[0], updated[k]))
                else:
                    f1.write(l)
    param['parameter_file']=new_parameter_file


################################################################################################################################################################################
################################################################################################################################################################################
####       Running Functions  ################################################################################################################################################
################################################################################################################################################################################
################################################################################################################################################################################
def addInputOutputFiles(param,input_files,output_files):
    #store information on which file type to work and what output file to use
    param['input_files'] = input_files
    param['output_files'] = output_files

    #check if the specified input files actually exist otherwise throw an error
    if not param.has_key(param['input_files']):
        writeLog('The input file type you specified does not exist')
        sys.exit(0)
        
    #check if the key for the specified output files actually exist otherwise create them
    if not param.has_key(param['output_files']) and param['output_files']!='':
         param[param['output_files']]=['']*param['num_samples']

def check_if_queue_finished_sucessful(param):
  # check how many ended successfully
    s_noend = []
    for idx in range(param['num_samples']):
        logfile = open(param['working_dir']+'results/log/'+param['stub'][idx]+'.log')
        lines_end = [line for line in logfile.readlines() if 'ENDING %s |' %(param['current_flag']) in line.rstrip()]
        #check if the job was successfully finished
        if len(lines_end) == 1:
            #check if there is suppose to be a output file and if that is the case store the location of the output file
            if  param['output_files']!='':
                retval=lines_end[0].split('|')[1].strip()
                if len(retval.split(';'))==1:
                    param[param['output_files']][idx]=retval
                else:
                    param[param['output_files']][idx]=retval.split(';')[0]
                    param[param['output_files']+'2'][idx]=retval.split(';')[1]    
            #add an entry into the run log that shows that the current job has been finished successfully
            param['run_log'][-1][idx]=True
        else:
            if  param['output_files']!='': param[param['output_files']][idx]=''
            s_noend.append(param['stub'][idx])
            
    sucessful=True
    #output error if not all samples finished successfully otherwise indicate that this step was completed sucessfully
    if len(s_noend) > 0:
        writeLog('error in samples %s' %(';'.join([s for s in s_noend])),param)
        writeLog('\n',param)
        sucessful=False
    else:        
        writeLog(param['current_flag']+' successful!\n\n',param)
    return sucessful

def dumpParameters(param):
    #dumps the parameter file into a JSON object
    with open(param['working_dir']+'results/parameters.json', 'w') as f:
        json.dump(param, f)

        
def is_module_finished(param):
    #check if the current module was already run successfully before 
    if param['clean_run']:
       finished = False
    else:
        handle=open(param['log_handle'])
        success = [line for line in handle.readlines() if '%s successful!' %(param['current_flag']) in line.rstrip()]
        handle.close()
        finished = len(success)==1
    return finished  

def submit_job(param,py_file,input_files,output_files='',cores='1-8',environment='modules_python2.6'):
    #store information on which files to use as input and which ones as output    
    addInputOutputFiles(param,input_files,output_files)

    #get the current flag for the log file
    param['current_dir']  = py_file.split('.')[0]
    param['current_flag'] = param['current_dir'] + '_' + input_files
        
    if is_module_finished(param):
        #add coment into main log file
        writeLog('Skipping '+param['current_flag']+' since it was already successfully finished.\n',param)
        #and fetch the current working files
        for idx in range(param['num_samples']):
            logfile = open(param['working_dir']+'results/log/'+param['stub'][idx]+'.log')
            #check if the current module actually produces output files that are used subsequently. (i.e. tophat or the trimmer as opposed to fastqc or bamqc)
            if  param['output_files']!='':
                #if there are output files fetch the values and add them to the proper list of files (bam_files, fastq_files, ... )
                lines_end = [line for line in logfile.readlines() if 'ENDING %s | ' %(param['current_flag']) in line.rstrip()]
                retval=lines_end[0].split('|')[1].strip()
                #check if there are one or 2 working files returned. e.g. the trimmer will return 2 values, while most modules return only 1
                if len(retval.split(';'))==1:
                    param[param['output_files']][idx]=retval
                else:
                    param[param['output_files']][idx]=retval.split(';')[0]
                    param[param['output_files']+'2'][idx]=retval.split(';')[1]    
        #add in the run log that all samples finished successfully
        param['run_log'].append([True]*param['num_samples'])
        param['run_log_headers'].append(param['current_flag'])
                
    else:
        writeLog('Starting '+param['current_flag']+'\n',param)
        #create current working directory
        param['module_dir']=param['working_dir']+'results/'+param['current_dir']+'/'
        if not os.path.exists(param['module_dir']):
            os.makedirs(param['module_dir'])

        #write all parameters to file so the single subnodes can use it
        dumpParameters(param)

        #generate a job id
        job_id='RNASeq_'+str(random.randint(1,1000))
 
        #submit all qsub jobs
        for index in range(param['num_samples']): 
            #check if the current sample was sucessfully processed before, otherwise skip the subsequent steps 
            if param['run_log'][-1][index]:
                if param['run_single_cpu']:
                    single_cpu_module.run_single_job(index, param, job_id, py_file, cores)    
                else:
                    qsub_module.submit_jobs(index, param, py_file, job_id, cores, environment)

        #wait for qsub scripts to finish
        if not param['run_single_cpu']:
            qsub_module.wait_for_qsub(param, job_id)

        #start a new log column for the current batch of jobs
        param['run_log'].append([False]*param['num_samples'])
        param['run_log_headers'].append(param['current_flag'])

        #check if all jobs finished successful
        check_if_queue_finished_sucessful(param)         
        
        writeLog('++++++++++++++++++++++++++++\n\n',param)

def writeLog(string,param):
    if param['verbose']: print string
    handle=open(param['log_handle'],'a')
    handle.write(string)
    handle.close()


################################################################################################################################################################################
################################################################################################################################################################################
####       Reporting Functions  ################################################################################################################################################
################################################################################################################################################################################
################################################################################################################################################################################

def copyFileAndLinkIt(param,paramKey,text):
    #copy file
    call='cp '+ param[paramKey]+' '+param['working_dir']+'report/'
    output,error = subprocess.Popen(call.split(),stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()
    #add file to html
    param['report'].write('<a href="'+param[paramKey].split('/')[-1]+'">'+text+'</a><br>')



def report_finish(param):
    param['report'].write('<style>table.fixed { table-layout:fixed; } table.fixed td { overflow: hidden; }#one-column-emphasis{font-family:"Lucida Sans Unicode", "Lucida Grande", Sans-Serif;font-size:12px;width:480px;text-align:left;border-collapse:collapse;margin:20px;}#one-column-emphasis th{font-size:14px;font-weight:normal;color:#039;padding:12px 15px;}#one-column-emphasis td{color:#669;border-top:1px solid #e8edff;padding:10px 15px;}.oce-first{background:#d0dafd;border-right:10px solid transparent;border-left:10px solid transparent;}#one-column-emphasis tr:hover td{color:#339;background:#eff2ff;}</style></body>\n')
    param['report'].close()

def report_run_log(param):
    param['report'].write('<center><br><h2>Running Log</h2>')   

    #create a table
    table=[]
    table.append([stub for stub in param['stub']])
    
    #links to the icons
    PASS='<img src="Icons/tick.png">'
    FAIL='<img src="Icons/error.png">'

    #create table    
    for idx in range(len(param['run_log_headers'])):
        line=[param['run_log_headers'][idx].title()]
        for i in range(param['num_samples']):
            status=FAIL
            if param['run_log'][idx][i]: status=PASS
            line=line+[status]
        table.append(line)
    
    #write the table as html     
    module_helper.writeHTMLtable(param,table,out=param['report'])

def report_start(param):
    writeLog('Starting to write the report ... \n',param)
    param['report']=open(param['working_dir']+'report/index.html','w')
    param['report'].write('<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd"><head><title></title></head><body>\n')    
    param['report'].write('<center><br><h1>RNASeq Report</h1>')
    copyFileAndLinkIt(param,'parameter_file','Parameter file')
    copyFileAndLinkIt(param,'raw_filenames','Raw files')
    
    #copy the pass/fail icons into the directory
    call='cp -R '+ param['scripts_dir']+'Icons '+param['working_dir']+'report/'
    output,error = subprocess.Popen(call.split(),stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()
 
    #report the run log in a table and show which samples passed/failed 
    report_run_log(param)    



