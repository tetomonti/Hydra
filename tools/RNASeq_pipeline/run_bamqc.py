"""
Usage: python2.7 bam_qc.py -i input_file.bam -o outdir
-h help
-o output_dir
-i input_file.bam				*[No default value]

"""


def extract_stats(input_file):
    #open bam file
    bam_file = pysam.Samfile(input_file, "rb")
    #counters
    total_aligned_reads=0; unique_aligned_reads=0; is_singleton=0; is_paired=0; is_proper_pair=0; is_unmapped=0 
    num_unique_mismatches=[0]*5; num_multiple_mismatches=[0.0]*5; num_multiread=[0.0]*20
    delet=False; insert=False; spliced=False
    reads_with_deletions=0; spliced_reads=0; reads_with_inserts=0; non_spliced_reads=0
    unique_reads_with_deletions=0; unique_spliced_reads=0; unique_reads_with_inserts=0; unique_non_spliced_reads=0

    #tag variables
    NH=0; NM=0
    idx=0
    for read in bam_file:
       #get all the relevant tags
        for tag in read.tags:
            if tag[0]=='NH': NH = tag[1]
            if tag[0]=='NM': NM = tag[1]
        #number of aligned reads 
        total_aligned_reads+=1
        unique_aligned_reads+=1/NH 
        
        #number of mismatches  
        if NH==1: 
            if NM>=4:
                num_unique_mismatches[4]=num_unique_mismatches[4]+1
            else:
                num_unique_mismatches[NM]=num_unique_mismatches[NM]+1
        else:
            if NM>=4:
                num_multiple_mismatches[4]=num_multiple_mismatches[4]+(1.0/float(NH))
            else:
                num_multiple_mismatches[NM]=num_multiple_mismatches[NM]+(1.0/float(NH))
        
        #number of multiple reads
        if NH>=20: 
           num_multiread[19]=num_multiread[19]+(1.0/float(NH))
        else:
           num_multiread[NH-1]=num_multiread[NH-1]+(1.0/float(NH))   
           
        #singletons, paired, proper paired, unmapped
        is_singleton   += int( not read.is_paired)
        is_paired      += int( read.is_paired)
        is_proper_pair += int( read.is_proper_pair)
        is_unmapped    += int( read.is_unmapped)
        
        #splicing, deletions, inserts
        spliced = 'N' in read.cigarstring
        insert  = 'I' in read.cigarstring
        delet   = 'D' in read.cigarstring
        
        #actual count
        spliced_reads       += int(spliced)
        spliced_reads       += int(spliced)
        non_spliced_reads   += int(not spliced)
        reads_with_deletions+= int(insert)
        reads_with_inserts  += int(delet)
             
        #counting reads that are aligned multiple times only once
        unique_spliced_reads       += int(spliced)/NH
        unique_non_spliced_reads   += int(not spliced)/NH
        unique_reads_with_deletions+= int(insert)/NH
        unique_reads_with_inserts  += int(delet)/NH 
        if idx % 1000000 == 0:
            print str(idx)+' reads done'
        idx+=1    
    bam_file.close()
    
    stats=dict()
    stats['total_aligned_reads']         = total_aligned_reads
    stats['unique_aligned_reads']        = unique_aligned_reads
    stats['is_singleton']                = is_singleton
    stats['is_paired']                   = is_paired
    stats['is_proper_pair']              = is_proper_pair
    stats['is_unmapped']                 = is_unmapped
    stats['num_unique_mismatches']       = num_unique_mismatches
    stats['num_multiple_mismatches']     = num_multiple_mismatches
    stats['num_multiread']               = num_multiread
    stats['spliced_reads']               = spliced_reads
    stats['non_spliced_reads']           = non_spliced_reads
    stats['reads_with_inserts']          = reads_with_inserts
    stats['reads_with_deletions']        = reads_with_deletions
    stats['unique_spliced_reads']        = unique_spliced_reads
    stats['unique_non_spliced_reads']    = unique_non_spliced_reads
    stats['unique_reads_with_inserts']   = unique_reads_with_inserts
    stats['unique_reads_with_deletions'] = unique_reads_with_deletions    
    return stats
    
    
def output_stats(stats, output_dir):
    #write all stats into a file
    handle=open(output_dir+'output.txt','w')
    handle.write('total_aligned_reads \t'+str(stats['total_aligned_reads'])+'\n')
    handle.write('unique_aligned_reads \t'+str(stats['unique_aligned_reads'])+'\n')
    handle.write('is_singleton \t'+str(stats['is_singleton'])+'\n')
    handle.write('is_paired \t'+str(stats['is_paired'])+'\n')
    handle.write('is_proper_pair \t'+str(stats['is_proper_pair'])+'\n')
    handle.write('is_unmapped \t'+str(stats['is_unmapped'])+'\n')
    for i in range(len(stats['num_unique_mismatches'])): handle.write('num_unique_mismatches '+str(i)+'\t'+str(stats['num_unique_mismatches'][i])+'\n')
    for i in range(len(stats['num_multiple_mismatches'])): handle.write('num_multiple_mismatches '+str(i)+'\t'+str(stats['num_multiple_mismatches'][i])+'\n')
    for i in range(len(stats['num_multiread'])): handle.write('num_multiread '+str(i+1)+'\t'+str(stats['num_multiread'][i])+'\n')
    handle.write('spliced_reads \t'+str(stats['spliced_reads'])+'\n')
    handle.write('non_spliced_reads \t'+str(stats['non_spliced_reads'])+'\n')
    handle.write('reads_with_inserts \t'+str(stats['reads_with_inserts'])+'\n')
    handle.write('reads_with_deletions \t'+str(stats['reads_with_deletions'])+'\n')
    handle.write('unique_spliced_reads \t'+str(stats['unique_spliced_reads'])+'\n')
    handle.write('unique_non_spliced_reads \t'+str(stats['unique_non_spliced_reads'])+'\n')
    handle.write('unique_reads_with_inserts \t'+str(stats['unique_reads_with_inserts'])+'\n')
    handle.write('unique_reads_with_deletions \t'+str(stats['unique_reads_with_deletions'])+'\n')
    handle.close()



def plot_mul_alignments(stats,output_dir):
    fig, ax = plt.subplots()
    index = np.arange(len(stats['num_multiread']))
    bar_width = 0.8
    opacity = 0.4
    val=[math.log(stat+1,10) for stat in stats['num_multiread']]
    rects1 = plt.bar(index, val, bar_width,
                     alpha=opacity,
                     color='b',
                     label='Number of alignements ')
    plt.xlabel('Number of alignments')
    plt.ylabel('Counts (log10)')
    plt.title('Distribution of reads with multiple alignments')
    ticks=[str(i+1)  for i in range(len(stats['num_multiread']))]
    ticks[len(ticks)-1]=ticks[len(ticks)-1]+'+'
    plt.xticks(index + bar_width, ticks)
    plt.tight_layout()
    pylab.savefig(output_dir+'multiple_alignments.png')


def plot_num_unique_mismatches(stats, output_dir):
    fig, ax = plt.subplots()
    index = np.arange(len(stats['num_unique_mismatches']))
    bar_width = 0.8
    opacity = 0.4
    val=[math.log(stat+1,10) for stat in stats['num_unique_mismatches']]
    rects1 = plt.bar(index, val, bar_width,
                     alpha=opacity,
                     color='b')
    plt.xlabel('Number of mismatches in uniquely aligned samples')
    plt.ylabel('Counts (log10)')
    plt.title('Distribution of mismatches in reads with unique alignments')
    ticks=[str(i)  for i in range(len(stats['num_unique_mismatches']))]
    ticks[len(ticks)-1]=ticks[len(ticks)-1]+'+'
    plt.xticks(index + bar_width, ticks)
    plt.tight_layout()
    pylab.savefig(output_dir+'num_unique_mismatches.png')

def num_multiple_mismatches(stats,output_dir):
    fig, ax = plt.subplots()
    index = np.arange(len(stats['num_multiple_mismatches']))
    bar_width = 0.8
    opacity = 0.4
    val=[math.log(stat+1,10) for stat in stats['num_multiple_mismatches']]
    rects1 = plt.bar(index, val, bar_width,
                     alpha=opacity,
                     color='b')
    plt.xlabel('Number of mismatches in multiple aligned samples')
    plt.ylabel('Counts (log10)')
    plt.title('Distribution of mismatches in reads with multiple alignments')
    ticks=[str(i)  for i in range(len(stats['num_multiple_mismatches']))]
    ticks[len(ticks)-1]=ticks[len(ticks)-1]+'+'
    plt.xticks(index + bar_width, ticks)
    plt.tight_layout()
    pylab.savefig(output_dir+'num_multiple_mismatches.png')



def create_html(stats,output_dir):
    handle=open(output_dir+'sample_stats.html','w')
    
    #output a table with all the counts
    handle.write('<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd"><head><title></title></head><body>\n')    
    handle.write('<center><br><h1>Sample overview</h1>')
    
    #table
    handle.write('<table id="one-column-emphasis">\n')
    handle.write('<thead><tr><th> </th><th>Count</th><th>Percentage</th></tr></thead>\n')
  
    #total number + unique / multiple aligned
    handle.write('<tr><td>Total number of aligned reads</td><td>'+str(int(stats['total_aligned_reads']))+'</td><td>'+str(100*round(float(stats['total_aligned_reads'])/float(stats['total_aligned_reads']),3))+'% </td></tr>\n')
    handle.write('<tr><td>Number of uniquely aligned reads</td><td>'+str(int(stats['num_multiread'][0]))+'</td><td>'+str(100*round(float(stats['num_multiread'][0])/float(stats['total_aligned_reads']),3))+'% </td></tr>\n')
    multi_read=stats['total_aligned_reads']-stats['num_multiread'][0]
    handle.write('<tr><td>Number of multiple aligned reads</td><td>'+str(int(multi_read))+'</td><td>'+str(100*round(float(multi_read)/float(stats['total_aligned_reads']),3))+'% </td></tr>\n')
    handle.write('<tr>  <td></td><td>  </td><td>  </td></tr>\n')
  
    #mismatches within uniquely aligned
    handle.write('<tr><td>Number of perfect matches within uniquely aligned reads</td><td>'+str(int(stats['num_unique_mismatches'][0]))+'</td><td>'+str(100*round(float(stats['num_unique_mismatches'][0])/float(stats['total_aligned_reads']),3))+'% </td></tr>\n')
    uniq_read_multi_mm=stats['num_multiread'][0]-stats['num_unique_mismatches'][0]
    handle.write('<tr><td>Number of uniquely aligned reads with mismatches</td><td>'+str(int(uniq_read_multi_mm))+'</td><td>'+str(100*round(float(uniq_read_multi_mm)/float(stats['total_aligned_reads']),3))+'% </td></tr>\n')
    handle.write('<tr>  <td></td><td>  </td>  <td>  </td></tr>\n')
    
    #mismatches within uniquely aligned
    handle.write('<tr><td>Number of perfect matches within multiple aligned reads</td><td>'+str(int(stats['num_multiple_mismatches'][0]))+'</td><td>'+str(100*round(float(stats['num_multiple_mismatches'][0])/float(stats['total_aligned_reads']),3))+'% </td></tr>\n')
    mul_read_multi_mm=multi_read-stats['num_multiple_mismatches'][0]
    handle.write('<tr><td>Number of multiple aligned reads with mismatches</td><td>'+str(int(mul_read_multi_mm))+'</td><td>'+str(100*round(float(mul_read_multi_mm)/float(stats['total_aligned_reads']),3))+'% </td></tr>\n')
    handle.write('<tr><td>   </td><td>   </td><td>   </td></tr>\n')

    #paired / singleton / ...
    handle.write('<tr><td>Number of singleton reads</td><td>'+str(stats['is_singleton'])+'</td><td>'+str(100*round(float(stats['is_singleton'])/float(stats['total_aligned_reads']),3))+'% </td></tr>\n')
    handle.write('<tr><td>Number of paired reads</td><td>'+str(stats['is_paired'])+'</td><td>'+str(100*round(float(stats['is_paired'])/float(stats['total_aligned_reads']),3))+'% </td></tr>\n')
    handle.write('<tr><td>Number of proper paired reads</td><td>'+str(stats['is_proper_pair'])+'</td><td>'+str(100*round(float(stats['is_proper_pair'])/float(stats['total_aligned_reads']),3))+'% </td></tr>\n')
    handle.write('<tr><td>Number of unmapped reads</td><td>'+str(stats['is_unmapped'])+'</td><td>'+str(100*round(float(stats['is_unmapped'])/float(stats['total_aligned_reads']),3))+'% </td></tr>\n')
    handle.write('<tr><td>   </td><td>   </td><td>   </td></tr>\n')
        
    #spliced / inserts / deletions
    handle.write('<tr><td>Number of spliced reads</td><td>'+str(stats['spliced_reads'])+'</td><td>'+str(100*round(float(stats['spliced_reads'])/float(stats['total_aligned_reads']),3))+'% </td></tr>\n')
    handle.write('<tr><td>Number of reads with inserts</td><td>'+str(stats['reads_with_inserts'])+'</td><td>'+str(100*round(float(stats['reads_with_inserts'])/float(stats['total_aligned_reads']),3))+'% </td></tr>\n')
    handle.write('<tr><td>Number of reads with deletions</td><td>'+str(stats['reads_with_deletions'])+'</td><td>'+str(100*round(float(stats['reads_with_deletions'])/float(stats['total_aligned_reads']),3))+'% </td></tr>\n')      
    handle.write('</table><br><br><br><br>\n')
    
    #add figures
    handle.write('<img src="multiple_alignments.png" alt="multiple_alignments"><br><br><br><br>\n')
    handle.write('<img src="num_unique_mismatches.png" alt="num_unique_mismatches"><br><br><br><br>\n')
    handle.write('<img src="num_multiple_mismatches.png" alt="num_multiple_mismatches"><center><br><br><br><br>\n\n\n')
    handle.write('<style>#one-column-emphasis{font-family:"Lucida Sans Unicode", "Lucida Grande", Sans-Serif;font-size:12px;width:480px;text-align:left;border-collapse:collapse;margin:20px;}#one-column-emphasis th{font-size:14px;font-weight:normal;color:#039;padding:12px 15px;}#one-column-emphasis td{color:#669;border-top:1px solid #e8edff;padding:10px 15px;}.oce-first{background:#d0dafd;border-right:10px solid transparent;border-left:10px solid transparent;}#one-column-emphasis tr:hover td{color:#339;background:#eff2ff;}</style></body>\n')

    handle.close()
    

def make_report(stats,output_dir):
    plot_mul_alignments(stats,output_dir)
    plot_num_unique_mismatches(stats,output_dir)
    num_multiple_mismatches(stats,output_dir)
    create_html(stats,output_dir)     
       
        
if __name__ == "__main__":
## Import modules
    import pysam, os, sys, getopt, json
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import math, pylab
    
## Check arguments
    if len(sys.argv)<5:
        print __doc__
        sys.exit(0)
    optlist,cmdlist = getopt.getopt(sys.argv[1:],'hi:o:') 

    for opt in optlist:
        if opt[0]=='-h': print __doc__; sys.exit(0)
        if opt[0]=='-i': input_file=opt[1]
        if opt[0]=='-o': output_dir=opt[1] 
        
    #extract stats from bam file
    stats=extract_stats(input_file)    
    
    #dump stats into a text file 
    output_stats(stats, output_dir)    
    
    #create a report for a single sample
    make_report(stats, output_dir)
    
    #dump stats into a json file
    with open(output_dir+'stats.json', 'w') as f:
        json.dump(stats, f)
    
    
    
