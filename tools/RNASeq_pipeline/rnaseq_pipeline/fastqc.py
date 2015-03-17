import os, subprocess, csv
import module_helper



def copyFiles(param, input_files):
    #if there is no fastqc directory in the report make one
    param['fastqc_dir']=param['working_dir']+'report/fastqc/'
    #Create a report fastqc subdirectory for each sample
    for s in param['stub']:
     if not os.path.exists(param['fastqc_dir']+'/'+s):
        os.makedirs(param['fastqc_dir']+'/'+s)
    
    #reconstruct file names
    param['fastqc_stub']=[stb+'/'+fn.split('/')[-1] for stb,fn in zip(param['stub'],param[input_files])]
    if param['paired']:
        param['fastqc_stub']=param['fastqc_stub']+[stb+'/'+fn.split('/')[-1] for stb,fn in zip(param['stub'],param[input_files+'2'])]
    
    #for fn in param['fastqc_stub']:
    # print fn

    #fix the suffix
    param['fastqc_stub']=[fn.replace('.fastq.gz','_fastqc') for fn in param['fastqc_stub']]
        
    #use only the files that exist
    fqc_dir=param['working_dir']+'results/fastqc/'
    #for fn in param['fastqc_stub']:
     #print fn
    param['fastqc_stub']=[fn for fn in param['fastqc_stub'] if os.path.exists(fqc_dir+fn)]
    #print param['fastqc_stub']
    #copy the unpacked directories   
    for fastqc_file in param['fastqc_stub']:
        call='cp -R '+ fqc_dir +fastqc_file+'/ '+param['fastqc_dir']+fastqc_file.split('/')[0]
#        print call
        output,error = subprocess.Popen(call.split(),stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()

    #cut off the suffix
    param['fastqc_stub']=[fn.replace('_fastqc','') for fn in param['fastqc_stub']]



def createOverviewTable(param):
    
#    print['Summary files']+['<a href="fastqc/'+stub+'_fastqc/fastqc_data.txt">raw</a>' for stub in param['fastqc_stub']]
#    print ['Full report']+['<a href="fastqc/'+stub+'_fastqc/fastqc_report.html"><img src="Icons/fastqc_icon.png"></a>' for stub in param['fastqc_stub']] 
    
    #create a table
    table=[]
    #put in headers
    table.append([stub for stub in param['fastqc_stub']])
    #link to summary files
    table.append(['Summary files']+['<a href="fastqc/'+stub+'_fastqc/fastqc_data.txt">raw</a>' for stub in param['fastqc_stub']])
    #link to overview files
    table.append(['Full report']+['<a href="fastqc/'+stub+'_fastqc/fastqc_report.html"><img src="Icons/fastqc_icon.png"></a>' for stub in param['fastqc_stub']])
    #extract check marks
    table=table+extractTables(param)    
    #write the table as html     
    module_helper.writeHTMLtable(param,table,out=param['report'])
     
def extractTables(param):
    #get the rownames
    csv_file = open(param['fastqc_dir']+param['fastqc_stub'][0]+'_fastqc/summary.txt')
    csv_reader = csv.reader(csv_file, delimiter='\t')
    #get the rownames    
    checks = [[row[1]] for row in csv_reader]
    csv_file.close()
  
    #links to the icons
    PASS='<img src="Icons/tick.png">'
    FAIL='<img src="Icons/error.png">'
    WARN='<img src="Icons/warning.png">'

    #get the values for each sample (icons for pass, faile or warning)
    for idx in range(len(param['fastqc_stub'])):
        csv_file = open(param['fastqc_dir']+param['fastqc_stub'][idx]+'_fastqc/summary.txt')
        overview_file='fastqc/'+param['fastqc_stub'][idx]+'_fastqc/fastqc_report.html#M'
        csv_reader = csv.reader(csv_file, delimiter='\t')

        i=0
        for row in csv_reader:
            cell='<a href="'+overview_file+str(i)+'">'
            if row[0]=='PASS': cell=cell+PASS
            if row[0]=='FAIL': cell=cell+FAIL
            if row[0]=='WARN': cell=cell+WARN
            cell=cell+'</a>'
            checks[i].append(cell)
            i+=1
        csv_file.close()
    return checks


def readRawfastqc(param,input_files):
    #print param['fastqc_dir']
    #print param['fastqc_stub']
    
    summary_files = [param['fastqc_dir']+param['fastqc_stub'][idx]+'_fastqc'+'/summary.txt' for idx in range(len(param['fastqc_stub']))]
    
    print "Summary files .."
    print summary_files

    fastqc=dict()
    #add entries into fastqc dictionary
    fh=open(summary_files[0])
    fastqc = dict([(name.split('\t')[1].strip(), []) for name in fh.readlines()])
    fh.close()

    #fill fastqc dictionary with information from the summary files
    for sum_file in summary_files:
        fh=open(sum_file)
        for name in (fh.readlines()):
            fastqc[name.split('\t')[1].rstrip()].append(name.split('\t')[0].strip())
        fh.close()
    
    key_list = fastqc.keys()
    
    #fill fastqc dictionary with information from the data file
    data_files = [param['fastqc_dir']+param['fastqc_stub'][idx]+'_fastqc'+'/fastqc_data.txt' for idx in range(len(param['fastqc_stub']))]
    labels = ['Encoding', 'Total Sequences', 'Filtered Sequences', 'Sequence length',	'%GC']
    #print "Data files .."
    #print data_files 
    fastqc.update(dict([(l, []) for l in labels]))
    key_list.extend(labels)
    
    for d_file in data_files:
        fh=open(d_file)
        for name in (fh.readlines()):
            if name.split('\t')[0].strip() in fastqc.keys():
                fastqc[name.split('\t')[0].strip()].append(name.split('\t')[1].rstrip())
        fh.close()
    
    param['fast_qc_summary']=fastqc

    # write overview file 
    fh=open(param['fastqc_dir']+input_files+'overview.txt','w')
    fh.write(' \t'+'\t'.join(param['fastqc_stub'])+'\n')
    for nam in key_list: fh.write(nam+'\t'+'\t'.join([str(vv) for vv in param['fast_qc_summary'][nam]])+'\n') 
    fh.close()



def plotNumberOfReads(param,input_files):
    import numpy as np
    import matplotlib.pyplot as plt
    import math, pylab
    #extract number of reads, these are needed not only here but also for the bamqc
    #for idx in range(len(param['fastqc_stub'])):
     #print idx
    summary_files = [param['fastqc_dir']+param['fastqc_stub'][idx]+'_fastqc/fastqc_data.txt' for idx in range(len(param['fastqc_stub']))]
    #print summary_files

    num_total_reads=[open(sum_file).readlines()[6].split('\t')[1].strip().rstrip() for sum_file in summary_files]
    num_total_reads=[int(num) for num in num_total_reads]
    
    #if we deal with unpaired data the total number of reads used downstream is equivalent to the fastqc reads, but for paired it is the sum of the two paired files
    #print range(param['num_samples'])
    #print num_total_reads
    #for idx in range(param['num_samples']):
    #print idx
    #print param[input_files][0]
    #print param[input_files][1]
#    print param['fastqc_stub']
    if not param['paired']:
        param['num_total_reads']=num_total_reads
    else:
        param['num_total_reads']=[0]*param['num_samples']
        #samples_1=[param[input_files][idx].split('/')[-1].replace('.fastq','').replace('.fq','').replace('.gz','') for idx in range(param['num_samples'])]
        #samples_2=[param[input_files][idx].split('/')[-1].replace('.fastq','').replace('.fq','').replace('.gz','') for idx in range(param['num_samples'])]
        samples_1=[param['fastqc_stub'][0]]
        samples_2=[param['fastqc_stub'][1]]
        
        #print samples_1
        #print samples_2
        #print param['stub']
#        print param['fastqc_stub'] 
#        print param['num_samples']
        #print num_total_reads 
        for idx in range(param['num_samples']):
        #    print idx
        #    index1=[i for i in range(len(param['fastqc_stub'])) if param['fastqc_stub'][i]==samples_1[idx]][0]
        #    index2=[i for i in range(len(param['fastqc_stub'])) if param['fastqc_stub'][i]==samples_2[idx]][0]           
            param['num_total_reads'][idx]=num_total_reads[idx]+num_total_reads[idx+param['num_samples']]
    #    print param['num_total_reads']
        
    #create plot 
    fig, ax = plt.subplots()
    fig.set_size_inches(3+len(param['fastqc_stub'])*0.4,8)
    index = np.arange(len(num_total_reads))
    bar_width = 0.8
    opacity = 0.4
    rects1 = plt.bar(index, num_total_reads, bar_width,
                     alpha=opacity,
                     color='b')
    plt.xlabel('Samples')
    plt.ylabel('Total number of reads')
    plt.title('Total number of reads across samples')
    ticks=param['fastqc_stub']
    plt.xticks(index + bar_width, ticks,rotation='vertical')
    plt.tight_layout()
    
    #put it into the report
    filename='report/fastqc/'+input_files+'total_reads.png'
    pylab.savefig(param['working_dir']+filename)
    param['report'].write('<img src="fastqc/'+input_files+'total_reads.png" alt="total number of reads"><br><br>\n')

def plotGCContent(param,input_files):
    import numpy as np
    import matplotlib.pyplot as plt
    import math, pylab

    #extract number of reads, these are needed not only here but also for the bamqc
    summary_files = [param['fastqc_dir']+param['fastqc_stub'][idx]+'_fastqc/fastqc_data.txt' for idx in range(len(param['fastqc_stub']))]
    gc_content=[open(sum_file).readlines()[9].split('\t')[1].strip().rstrip() for sum_file in summary_files]
    gc_content=[int(num) for num in gc_content]
    
    #create plot 
    fig, ax = plt.subplots()
    fig.set_size_inches(3+len(param['fastqc_stub'])*0.4,8)
    index = np.arange(len(gc_content))
    bar_width = 0.8
    opacity = 0.4
    rects1 = plt.bar(index, gc_content, bar_width,
                     alpha=opacity,
                     color='b')
    plt.xlabel('Samples')
    plt.ylabel('%GC content')
    plt.title('GC content across samples')
    ticks=param['fastqc_stub']
    plt.xticks(index + bar_width, ticks,rotation='vertical')
    plt.tight_layout()
    ax.set_ylim(0, 100)
    
    #put it into the report
    filename='report/fastqc/'+input_files+'gc_content.png'
    pylab.savefig(param['working_dir']+filename)
    param['report'].write('<img src="fastqc/'+input_files+'gc_content.png" alt="GC content"><br><br>\n')


def plotAgvLengthOfReads(param):
    import numpy as np
    import matplotlib.pyplot as plt
    import math, pylab

    #extract number of reads, these are needed not only here but also for the bamqc
    summary_files = [param['fastqc_dir']+param['fastqc_stub'][idx]+'_fastqc/fastqc_data.txt' for idx in range(len(param['fastqc_stub']))]
    seq_length=[open(sum_file).readlines()[8].split('\t')[1].strip().rstrip() for sum_file in summary_files]
    seq_length=[int(num) for num in seq_length]
    
    #create plot 
    fig, ax = plt.subplots()
    fig.set_size_inches(3+len(param['fastqc_stub'])*0.4,6)
    index = np.arange(len(seq_length))
    bar_width = 0.8
    opacity = 0.4
    rects1 = plt.bar(index, seq_length, bar_width,
                     alpha=opacity,
                     color='b')
    plt.xlabel('Samples')
    plt.ylabel('Read length')
    plt.title('Read length across samples')
    ticks=param['fastqc_stub']
    plt.xticks(index + bar_width, ticks,rotation='vertical')
    plt.tight_layout()
    
    #put it into the report
    filename='report/fastqc/read_length.png'
    pylab.savefig(param['working_dir']+filename)
    param['report'].write('<img src="fastqc/read_length.png" alt="read length"><br><br>\n')


def report(param,input_files='fastq_files',header='FastQC results'):  
    #assemble the full fastqc report
    copyFiles(param,input_files)
    param['num_total_reads']=[0]*param['num_samples']
    if len(param['fastqc_stub'])>0:
        param['report'].write('<center><br><h2>'+header+'</h2>')    
        readRawfastqc(param,input_files)
        createOverviewTable(param)
        param['report'].write('<a href="fastqc/'+input_files+'overview.txt">Table as tab delimited file</a><br><br><br>')
        plotNumberOfReads(param,input_files)
        plotGCContent(param,input_files)
        #plotAgvLengthOfReads(param)
    else:
        param['report'].write('There were no results to show.')    


def init(param):
   module_helper.checkParameter(param,key='fastqc_exec',dType=str,checkFile=True)

def run_fastqc(filename,param):
    out_dir = param['module_dir']+'/'+param['stub'][param['file_index']]
    if not os.path.exists(out_dir):
       os.makedirs(out_dir)
    call = param['fastqc_exec']+" '"+param[filename]+"' -o "+out_dir+' --extract'
    param['file_handle'].write('CALL: '+call+'\n')
    output,error = subprocess.Popen(call ,stdout = subprocess.PIPE,shell=True, stderr= subprocess.PIPE).communicate()
    param['file_handle'].write(error)
    param['file_handle'].write(output)
    if ('Analysis complete for' not in output):
       sys.exit()


if __name__ == "__main__":
    import subprocess, sys
    param=module_helper.initialize_module()

    #run fastqc
    run_fastqc('working_file',param)
    
    if not param['paired']:
        module_helper.wrapup_module(param)
    #calling it on the second fastq file if it is paired    
    else:
        run_fastqc('working_file2',param)
        module_helper.wrapup_module(param) 


