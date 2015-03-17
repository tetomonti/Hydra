import module_helper, helper, os, subprocess

def init(param):
   module_helper.checkParameter(param,key='cufflinks_exec',dType=str)
   module_helper.checkParameter(param,key='cufflinks_G',dType=str,checkFile=True)
   module_helper.checkParameter(param,key='cufflinks_b',dType=str,checkFile=True)
   module_helper.checkParameter(param,key='cufflinks_library_type',dType=str)
   module_helper.checkParameter(param,key='cufflinks_compatible_hits',dType=str)
   module_helper.checkParameter(param,key='cufflinks_total_hits',dType=str)
   module_helper.checkParameter(param,key='cufflinks_N',dType=str)
   module_helper.checkParameter(param,key='cufflinks_u',dType=str)   
   
   

def report(param):
    param['report'].write('<center><br><br><h2>Cufflinks results</h2>')    

    #make cufflinks directory
    cufflinks_dir=param['working_dir']+'report/cufflinks/'
    if not os.path.exists(cufflinks_dir):
        os.makedirs(cufflinks_dir)

    #output_phenotype_file
    module_helper.output_sample_info(param)
        
    #run R script that creates a PCA
    counts=param['working_dir']+'deliverables/cufflinks_counts_fpkm.txt'
    cmd=[param['Rscript_exec'],'cufflinks_pca.R','-c',counts,'-a',param['pheno_file'],'-o',param['working_dir']]
    output,error = subprocess.Popen(cmd ,stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()
    param['report'].write('<a href="cufflinks/pca.html">PCA</a>')
        
    
def finalize(param,input_files='count_files'):
    helper.writeLog('Collecting Cufflinks FPKM data ... \n',param)
    #extracts the counts from the cufflinks output
    import csv
 
    #check which of these files are actually available
    working_files = [iFile for iFile in param[input_files] if iFile !='']
    
    if len(working_files)>0:
        #get gene annotation
        csv_file = open(working_files[0])
        csv_reader = csv.reader(csv_file, delimiter='\t')
        counts = [row[0]+'\t'+row[4] for row in csv_reader]
        csv_file.close()
    
        #get all the expression values
        header='ENS_ID\tSymbol'
        for idx in range(param['num_samples']):
            if param[input_files]!='':            
                header=header+'\t'+param['stub'][idx]
                csv_file = open(param[input_files][idx])
                csv_reader = csv.reader(csv_file, delimiter='\t')
                i=0
                for row in csv_reader:
                    counts[i] = counts[i]+'\t'+row[9]
                    i+=1 
                csv_file.close()
    
        #output the file
        out_handle=open(param['working_dir']+'deliverables/cufflinks_counts_fpkm.txt','w')
        out_handle.write(header+'\n')
        for i in range(len(counts)): out_handle.write(counts[i]+'\n')
        out_handle.close()
    else:
        print('Cufflinks was not run successfully on any of the files..\n')
    
    
if __name__ == "__main__":
    import subprocess, sys, os
    param=module_helper.initialize_module()

    #run create output directory
    outdir = param['module_dir']+param['stub'][param['file_index']]+'/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    #build cufflinks call:
    call = param['cufflinks_exec']+' '
    
    #optional parameters
    if param['cufflinks_compatible_hits']=='active':
        call=call+'--compatible-hits-norm '
    if param['cufflinks_N']=='active':
        call=call+'-N '
    if param['cufflinks_u']=='active':
        call=call+'-u '
    if param['cufflinks_total_hits']=='active':
        call=call+'--total-hits-norm '
        
    #finish call
    call=call+' --no-update-check --library-type '+param['cufflinks_library_type']+' -o '+outdir+' -p '+param['num_processors']+' -G '+param['cufflinks_G']+' '+ param['working_file']
    
    param['file_handle'].write('CALL: '+call+'\n')
    output,error = subprocess.Popen(call.split() ,stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()
    param['file_handle'].write(output)
    
    #printing the error blows the log file up to several megabytes
    #param['file_handle'].write(error)
     
    #error handling
    if not os.path.exists(outdir+'genes.fpkm_tracking'):
        param['file_handle'].write('Error there was no cuffinks output\n') 
        sys.exit(0)
    
    #wrap up and return the current workingfile
    module_helper.wrapup_module(param,[outdir+'genes.fpkm_tracking']) 


