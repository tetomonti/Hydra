from rnaseq_pipeline import module_helper, helper
import subprocess, os

def init(param):
   module_helper.checkParameter(param,key='featureCount_exec',dType=str)
   #module_helper.checkParameter(param,key='sam_exec',dType=str)
   module_helper.checkParameter(param,key='featureCount_s',dType=str)
   module_helper.checkParameter(param,key='featureCount_t',dType=str)
   #module_helper.checkParameter(param,key='featureCount_m',dType=str,allowed=['union', 'intersection-strict','intersection-nonempty'])
   module_helper.checkParameter(param,key='featureCount_id',dType=str)
   module_helper.checkParameter(param,key='featureCount_gft',dType=str,checkFile=True)
   module_helper.checkParameter(param,key='featureCount_by_meta',dType=bool)
   module_helper.checkParameter(param,key='Rscript_exec',dType=str)
   
def createESet(count_file,pheno_file,param):
    #create a Bioconductor ExpresionSet
    cmd=[param['Rscript_exec'],'createRawCountESet.R','-c',count_file,'-a',pheno_file,'-o',param['working_dir']]
    output,error = subprocess.Popen(cmd ,stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()
    helper.writeLog('Creating ESet ... \n',param)
    

def copyAndReadStatFiles(param):
    #if there is no featureCount directory in the report make one
    featureCount_dir=param['working_dir']+'report/featureCount/'
    if not os.path.exists(featureCount_dir):
        os.makedirs(featureCount_dir)
 
    #get the files that are actually in the output directory       
    call=['cp','-R',param['working_dir']+'results/featureCount/',param['working_dir']+'report/']
    output,error = subprocess.Popen(call,stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()
    print error
    print output
    
    featureCount_file=param['working_dir']+'results/featureCount/featureCount_stats.txt'    
    #extract table
    table=[]
    fh=open(featureCount_file)
    #header
    table.append(fh.readlines()[0].rstrip().split('\t'))    
    table[0]=table[0][1:]
    fh.close()
   
    #total number of aligned reads
    tot_reads=['Total number of aligned reads']+param['bam_qc']['total_aligned_reads']
    table.append(tot_reads)    

    fh=open(featureCount_file)
    for line in (fh.readlines()[1:]):
       cur_line=line.rstrip().split('\t')
       perc=[cur_line[0]]+module_helper.getPercentage(cur_line[1:],tot_reads[1:],len(cur_line)-1)   
       table.append(perc)
              
    fh.close()
    return table
 


def report(param):
    #report only if there were actually results
    if os.path.exists(param['working_dir']+'deliverables/featureCount_raw_counts.txt'):
        param['report'].write('<center><br><br><br><br><h2>featureCount statistics</h2>')    
        table=copyAndReadStatFiles(param)
        module_helper.writeHTMLtable(param,table,out=param['report'],cell_width=80,fCol_width=150,deg=315)
        param['report'].write('<a href="featureCount/featureCount_stats.txt">featureCount statistics as tab delimited txt file</a>')
        


def finalize(param,input_files='count_files'):
    helper.writeLog('Collecting featureCount raw counts ... \n',param)
    #extracts the counts from the featureCount output
    import csv
 
    #check which of these files are actually available
    working_files = [iFile for iFile in param[input_files] if iFile !='']
    
    if len(working_files)>0:
        #get feature ID using the first column in the first file in the list of working files
        csv_file = open(working_files[0])
        csv_reader = csv.reader(csv_file, delimiter='\t')
 
        #For featureCount output, we want to skip the first two lines as they include the featureCount call and the headers which we don't want
        next(csv_reader,None)
        next(csv_reader,None)

        #Now start by taking the list of identifier, which is the first column in the file
        counts = [row[0] for row in csv_reader] 
        csv_file.close()
    
        #get all the expression values
        header='ID'
        for idx in range(param['num_samples']):
            if param[input_files]!='':            
                header=header+'\t'+param['stub'][idx]
                csv_file = open(param[input_files][idx])
                csv_reader = csv.reader(csv_file, delimiter='\t')
          
                #Here too we want to skip the first two lines, before getting the counts
                next(csv_reader,None)
                next(csv_reader,None)
                #Now start getting the counts (row[6]) and add in the ID (counts[i]) before it
                i=0
                for row in csv_reader:
                    counts[i] = counts[i]+'\t'+row[6]
                    i+=1 
                csv_file.close()
    
        #output the file
        out_file=param['working_dir']+'deliverables/featureCount_raw_counts.txt'
        out_handle=open(out_file,'w')
        out_handle.write(header+'\n')
        
        for i in range(len(counts)): out_handle.write(counts[i]+'\n')
        out_handle.close()
        
        #output_phenotype_file
        helper.writeLog('Writing phenotype data ... \n',param)
        module_helper.output_sample_info(param)
        
        #create an eSet:
        createESet(out_file,param['pheno_file'],param)
        
        #write summary stats        # featureCount does this on its own so we can just fetch each summary file
        #check which of these files are actually available
        working_files = [iFile+'.summary' for iFile in param[input_files] if iFile !='']
      
        if len(working_files)>0:
        #get Status column from summary file using the first column in the first file in the list of working files
         csv_file = open(working_files[0])
         csv_reader = csv.reader(csv_file, delimiter='\t')
        #Here, we want to skip the first line, as it simply points to the alignment file used when running featureCount
         next(csv_reader,None)
    
        #Now start by taking the list of identifier, which is the first column in the file
         entry = [row[0] for row in csv_reader]
         csv_file.close()

        #get all the summary stats for each sample
         header='Status'
      
         for idx in range(param['num_samples']):
            if param[input_files]!='':
              header=header+'\t'+param['stub'][idx]
             #Fetch the corresponding sample's summary file
              csv_file = open(param[input_files][idx]+'.summary')
              csv_reader = csv.reader(csv_file,delimiter='\t')
             #Again, we want to skip the first line, as it simply points to the alignment file used when running featureCount
              next(csv_reader,None)

             #Now start getting the stats (row[1]) and add in the Status (counts[i]) before it
              i=0
              for row in csv_reader:
                entry[i] = entry[i]+'\t'+row[1]
                i+=1
              csv_file.close()


        #output the file
         out_handle=open(param['working_dir']+'results/featureCount/featureCount_stats.txt','w')
         out_handle.write(header+'\n')
       
         for i in range(len(entry)): out_handle.write(entry[i]+'\n')
         out_handle.close()
        
    else:
        print('featureCount was not run successfully on any of the files..\n')
        
        
        
def main():
    import subprocess, sys, os
    param=module_helper.initialize_module()

    #build featureCount-count call:
    #If single-end then we don't need to add additional flags for paired-end data
    outfile_prefix=param['module_dir']+param['stub'][param['file_index']]
    if param['paired']:
      if param['featureCount_by_meta']:
        outfile = outfile_prefix+'_'+param['featureCount_id']+'.txt'
	call = [param['featureCount_exec'],'-P','-p','-t',param['featureCount_t'],'-s',param['featureCount_s'],'-g',param['featureCount_id'],'-a',param['featureCount_gft'],'-o',outfile,param['working_file']]
      else:
        outfile = outfile_prefix+'_'+param['featureCount_t']+'.txt'
        call = [param['featureCount_exec'],'-P','-p','-f','-t',param['featureCount_t'],'-s',param['featureCount_s'],'-g',param['featureCount_id'],'-a',param['featureCount_gft'],'-o',outfile,param['working_file']]
    else:
      if param['featureCount_by_meta']:
        outfile = outfile_prefix+'_'+param['featureCount_id']+'.txt'
	call = [param['featureCount_exec'],'-t',param['featureCount_t'],'-s',param['featureCount_s'],'-g',param['featureCount_id'],'-a',param['featureCount_gft'],'-o',outfile,param['working_file']]
      else:
	outfile = outfile_prefix+'_'+param['featureCount_t']+'.txt'
	call = [param['featureCount_exec'],'-f','-t',param['featureCount_t'],'-s',param['featureCount_s'],'-g',param['featureCount_id'],'-a',param['featureCount_gft'],'-o',outfile,param['working_file']]
    
    #print call
    #function calls
    param['file_handle'].write('CALL: '+' '.join(call)+'\n')
    output,error = subprocess.Popen(call ,stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()

    param['file_handle'].write(error)
    param['file_handle'].write(output)
        
    #error handling
    if not os.path.exists(outfile):
        param['file_handle'].write('featureCount run failed \n') 
        sys.exit(0)
    
    #write output
    #outfile=param['module_dir']+param['stub'][param['file_index']]+'.txt'
    #handle=open(outfile,'w')
    #handle.write(output)
    #handle.close()
    
    #wrap up and return the current workingfile
    module_helper.wrapup_module(param,[outfile]) 
     


