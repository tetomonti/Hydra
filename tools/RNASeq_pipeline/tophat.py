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
import module_helper

def init(param):
   module_helper.checkParameter(param,key='tophat_exec',dType=str)
   module_helper.checkParameter(param,key='tophat_index',dType=str)
   module_helper.checkParameter(param,key='tophat_qual',dType=str,optional=True)
   module_helper.checkParameter(param,key='tophat_N',dType=str)
   module_helper.checkParameter(param,key='tophat_gap_length',dType=str)
   module_helper.checkParameter(param,key='tophat_edit_dist',dType=str)
   module_helper.checkParameter(param,key='mate_inner_dist',dType=str)
   module_helper.checkParameter(param,key='mate_std_dev',dType=str)
   
   
     
if __name__ == "__main__":
    import subprocess, sys, os
    param=module_helper.initialize_module()

    #run create output directory
    outdir = param['module_dir']+param['stub'][param['file_index']]+'/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    #build tophat call:
    call = param['tophat_exec']+' '
    
    #resume option 
    if param['resume_module']:
        call=call+'-R '+outdir
    
    else:  
        #optional quality value
        if param['tophat_qual']!='none':
            call=call+param['tophat_qual']+' '
    
        #add paired parameters if paired
        if param['paired']:
           call = call + ' --mate-inner-dist ' + param['mate_inner_dist'] + ' --mate-std-dev ' + param['mate_std_dev'] + ' '
    
        call = call+'-o '+outdir+' -p '+param['num_processors']+' -N '+param['tophat_N']+' --read-gap-length '+param['tophat_gap_length']+' --read-edit-dist '+param['tophat_edit_dist']+' '+param['tophat_index']+' '+ param['working_file']

        
        #if paired add second working file
        if param['paired']:
            call = call+' '+param['working_file2']

       
    param['file_handle'].write('CALL: '+call+'\n')
    output,error = subprocess.Popen(call.split() ,stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()
    param['file_handle'].write(error)
    param['file_handle'].write(output)
       
    #error handling - checking if tophat logfile exists and if it does, check if tophat was run successful
    tophat_logfile = outdir+'logs/tophat.log'

    if not os.path.exists(tophat_logfile):
        sys.exit(0)
    else:
        log = open(tophat_logfile)
        lines_end = [line for line in log.readlines() if ' Run complete: ' in line.rstrip()]
        log.close()
        if len(lines_end) == 0:
            sys.exit(0)
    
    #wrap up and return the current workingfile
    module_helper.wrapup_module(param,[outdir+'accepted_hits.bam']) 




