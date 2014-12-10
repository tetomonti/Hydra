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

