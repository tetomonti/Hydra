
============================
Modification of the pipeline
============================

Overview
========

In short you need to add one file, which contains the code to run your new module on the command line 
and you need to change 2 additional files, one of them being the parameter file to include the module 
specific parameters and pipeline.py, which contains the workflow.

Detailed walkthrough
====================

1) Adding star specific parameters to the parameter file
--------------------------------------------------------

All parameters are contained in the parameter.txt file, which is a regular txt flat file. 
'#' are used to comment, so everything after '#' will be ignored by the pipeline. 
You can and are encouraged to add helpful comments and default parameters. 

.. code:: bash

   key = value  #comment

Parameters are defined as key and value pairs, which are separated by '='. The pipeline does not recognize 
white space so ' ' and tabs are removed, which mean one cannot define values that consist of words or values separated by spaces. 
Everything on the left of '=' is defined as key and everything on the right as value. 

The two most important parameters to add are the executable for Star and the pre-built index that corresponds to the proper species. 
Since Star is already included in our conda enviroment all we add is STAR as ``star_exec``. If you would like to use another version 
you can specify the directory and executable, but we do not offer support in this case (and therefore haven't tested it extensively), 
so you have to make sure it plays nicely with our current python2.7.9 version. Alternatively, you can read our documentation on how to 
deploy your own version of a program into the conda environment (see Developer's guide)

The second parameter ``star_index`` is a pre-built index, which has to be built only once, for which we provide the code in the supplement 
of this document. Star also allows the specification of a large variety of parameters. We are not changing these, but provide them in 
the parameter file in case someone would like to change it. The most important one here is ``outputSAMtype`` which we specify as 
``BAM_SortedByCoordinate`` by default Star outputs the aligned reads unsorted in the sam format. However, many tools require sorted reads 
and the next step in our pipeline requires .bam files, so it is more convenient to output them directly as sorted bam files.

.. code:: bash

   #######################################################################################
   STAR ALIGNER
   ######################################################################################
   run_star                  =   TRUE
   star_exec                 =   STAR     
   #The Star index index must be built prior to running star
   star_index                =   /restricted/projectnb/montilab-p/CBMrepositoryData/annot
   /hg19_ensembl_GRCh37_tophat/Ensembl/GRCh37/Sequence/star_index/

   #optional parameters (do not change unless you've read the documentation):
   outFilterType         = BySJout #default is BySJout
   outFilterMultimapNmax = 20    #max number of multiple alignments allowed (def: 20)
   alignSJoverhangMin    = 8     #minimum overhang for unannotated junctions (def: 8)
   alignSJDBoverhangMin  = 1     #minimum overhang for annotated junctions (def: 1)
   outFilterMismatchNmax = 999   #max number of mismatches per pair (def: 999, inactive)
   outFilterMismatchNoverLmax = 0.04  #max number of mismatches per pair relative to the 
   #read length(def: 0.04, 8 for paired reads with length 100)
   alignIntronMin        = 20  #minimum intron length
   alignIntronMax        = 1000000  #maximum intron length
   alignMatesGapMax      = 1000000  #maximum genomic distance between mates
   outputSAMtype         = BAM_SortedByCoordinate  #Can be BAM_SortedByCoordinate or BAM_unsorted

   
2) Creating the actual module file
----------------------------------

The most straight-forward way is to copy a module with similar functionality and modify it. In this case we could start from tophat.py, but here 
we start from scratch and go through step by step. First step is to load the helper libraries from the pipeline so we don't have to worry about 
parallelization and such:

.. code:: bash

	import rnaseq_pipeline.module_helper
	MODULE_HELPER = rnaseq_pipeline.module_helper

**Parameter initialization**

Next step is the initialization, where we need to initialize all parameters that we specified in the parameter.txt file. The pipeline has a 
helper function that allows to check if a key value pair is set and allows to cast to different data types.:

.. code:: bash

	check_parameter(param, key, dtype, allowed=[], checkfile=False, optional=False)   
	
``param`` ... is the parameter object that is provided by the pipeline, which acts as dynamic storage

``key`` ... key from the parameter file

``dtype`` ... data type for the value, which can be any python data type (str, int, float, ...). Of not is that 

``allowed`` ... is a list of allowed parameters, which is check against the value in the parameter file

``checkfile`` ... checks if the file physically exists, useful for annotations and the like

``optional`` ... in some cases this parameters are not needed mandatory and if not specified will be specified as blank - ''

We can use the ``check_parameter`` in the ``init()``function to check all the parameters. Of note here is that we are using all 
values as strings rather than numbers or integers, since we do not use this values for calculations, but only to build function calls later on:

.. code:: bash

	def init(param):
    		MODULE_HELPER.check_parameter(param, key='star_exec', dtype=str)
    		MODULE_HELPER.check_parameter(param, key='star_index', dtype=str, checkfile=True)
    		MODULE_HELPER.check_parameter(param, key='outFilterType', dtype=str)
    		MODULE_HELPER.check_parameter(param, key='outFilterMultimapNmax', dtype=str)
    		MODULE_HELPER.check_parameter(param, key='alignSJoverhangMin', dtype=str)
    		MODULE_HELPER.check_parameter(param, key='alignSJDBoverhangMin', dtype=str)
    		MODULE_HELPER.check_parameter(param, key='outFilterMismatchNmax', dtype=str)
    		MODULE_HELPER.check_parameter(param, key='outFilterMismatchNoverLmax', dtype=str)
    		MODULE_HELPER.check_parameter(param, key='alignIntronMin', dtype=str)
    		MODULE_HELPER.check_parameter(param, key='alignIntronMax', dtype=str)
    		MODULE_HELPER.check_parameter(param, key='alignMatesGapMax', dtype=str)
    		MODULE_HELPER.check_parameter(param, key='outputSAMtype', 
                                  allowed=['BAM_SortedByCoordinate',
                                           'BAM_unsorted'], 
                                  dtype=str)


**Main function to run on a single sample:**

With all the parameters initialized we can now write a function that builds a command-line call 
that is run on each sample. In the background, the pipeline handles all the scheduling and submission to the cluster, 
so all you have to worry about is building the actual function call. The entire script is run on each sample so we need 
to create a main function. (you can still define additional functions in the script, in this case we just put everything into the main:

.. code:: bash

	def main():
	
Before building the call we need to get all require variables, get the pointers to the files we want to work on and so on. 
There is a function that takes care of that: ``MODULE_HELPER.initialize_module()`` which should be called at the beginning of 
the main function.
	
.. code:: bash

	param = MODULE_HELPER.initialize_module()
	
There are several variables already initialized that make your life easier these include:

``param['module_dir']``    	-  full path to the current working directory, all output goes in here

``param['file_index']``    	- The index of the sample that is currently processed

``param['stub']``   	- An array of all output stubs as specify in the raw filename file

``param['outstub']`` 	- Output stub of the current working file, use that as part of your output

``param['working_file']`` 	- Current working file, on which the current tool should be run

``param['working_file2']``	- Current working file 2, for modules run on paired end seq runs, up to the alignment module

``param['file_handle']``	- Log file handle, which writes into the samples specific log file.

``param['paired']``	- Boolean flag that indicates whether this is a paired end seq run

Here we make a directory for each sample. Many tools have fixed output filenames. Having a directory for each sample avoids 
overwriting during parallelization

.. code:: bash

    #run create output directory
    outdir = param['module_dir']+param['outstub']+'/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

Next we build the command line call, which is basically just a list of all commands that are put into the command line. 
There are no spaces allowed, so the command line command is broken up into its pieces:

.. code:: bash

    call = [param['star_exec']]

    #add the directory where we have built the star index
    call.append('--genomeDir')
    call.append(param['star_index'])

    #add the number of processors to use
    call.append('runThreadN')
    call.append(param['num_processors'])

    #add all the optional parameters
    call.append('--outFilterType')
    call.append(param['outFilterType'])
    call.append('--outFilterMultimapNmax')
    call.append(param['outFilterMultimapNmax'])
    call.append('--alignSJoverhangMin')
    call.append(param['alignSJoverhangMin'])
    call.append('--alignSJDBoverhangMin')
    call.append(param['alignSJDBoverhangMin'])
    call.append('--outFilterMismatchNmax')
    call.append(param['outFilterMismatchNmax'])
    call.append('--outFilterMismatchNoverLmax')
    call.append(param['outFilterMismatchNoverLmax'])
    call.append('--alignIntronMin')
    call.append(param['alignIntronMin'])
    call.append('--alignIntronMax')
    call.append(param['alignIntronMax'])
    call.append('--alignMatesGapMax')
    call.append(param['alignMatesGapMax'])
    
    
We need to specify the output file type, most of the time this is going to be bam sorted by coordinate. 
Star changes the output file type 

.. code:: bash

	if (param['outputSAMtype'] == 'BAM_SortedByCoordinate'):
	        call.append('--outSAMtype')
	        call.append('BAM')
	        call.append('SortedByCoordinate')
	        outfile = 'Aligned.sortedByCoord.out.bam'
	elif (param['outputSAMtype'] == 'BAM_unsorted'):
		call.append('--outSAMtype')
        	call.append('BAM')
        	call.append('Unsorted')
        	outfile = 'Aligned.out.bam'
    	else:
    		outfile = 'Aligned.out.sam'
    	#add the proper output directories
    	call.append('--outFileNamePrefix')
    	call.append(outdir)


Using the list format to build the command allows for building alternative commands based on flags:

.. code:: bash

    #specify whether the fastq files are zipped
    call.append('--readFilesCommand')
    if param['zipped_fastq']:
        call.append('gunzip')
        call.append('-c')
    else:
        call.append('UncompressionCommand')


At the end we add the working files:

.. code:: bash

    #adding the files we want to work on
    call.append('--readFilesIn')
    call.append(param['working_file'])


The ``param['paired']`` flag indicates whether this is a paired end run. If it is we need to input the second working file as well.

.. code:: bash

    #if paired add second working file
    if param['paired']:
        call.append(param['working_file2'])

We found that it is good practice to output the complete function call into the log file. 
If there are errors in building the function call this lets you copy and paste it into the command line and find the bug much quicker. 
The function call is provided as list so we need to link that list using spaces:

.. code:: bash

    param['file_handle'].write('CALL: '+' '.join(call)+'\n')


We use subprocess to run the actual function call. This can be more sophisticated and include piping. 
For a more sophisticated example of this look into htseq.py.

.. code:: bash

    output, error = subprocess.Popen(call,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE).communicate()
                                     

Writing the output resulting from subprocess into the log file, so we have one place to look into. 
Some tools have a lot of output and it is better to leave them out. 

.. code:: bash

    param['file_handle'].write(error)
    param['file_handle'].write(output)


Finally, we need to wrap up the run, close the log file and so forth. 
Again there is already a function that does that for you. All you need to provide is the location of the output file(s). 
In this case we only provide one file since both pairs are already aligned after running Star, however, for the steps run 
before Star, such as cutadapt, you need to specify both working files:

.. code:: bash

    MODULE_HELPER.wrapup_module(param, [outdir+'OUTPUTFILENAME'])


**Reporting function:**

The star aligner does not have output that we would like to put into an html report, but other functions such as fastqc do. 
Such modules need to include a report function, which outputs into an html report. The structure is not quite as standardized as the initialization or running functions, but examples can be found in the fastqc.py and bamqc.py.

**Finalize function:**

Another function that is used only in a few function is a merge function that collects all the single results and creates a final matrix with all samples. htseq.py, cufflinks.py and featureCounts.py are examples for that.

3) Modifying setup.py
---------------------

Once we have created our star.py wrapper, we need to make sure that it gets installed with the RNASeq pipeline. 
For that we need to add the entry point into the setup.py file:
    
.. code:: bash

	entry_points={
        'console_scripts': [
            'run_bamqc=rnaseq_pipeline.bamqc:main',
            'run_cufflinks=rnaseq_pipeline.cufflinks:main',
            'run_cutadapt=rnaseq_pipeline.cutadapt:main',
            'run_fastqc=rnaseq_pipeline.fastqc:main',
            'run_featureCount=rnaseq_pipeline.featureCount:main',
            'run_htseq=rnaseq_pipeline.htseq:main',
            'run_matched_pairs=rnaseq_pipeline.matched_pairs:main',
            'run_tophat=rnaseq_pipeline.tophat:main',
            'run_star=rnaseq_pipeline.star:main']
    }
    
    
This provides us with the means to run 'run_star' in the command line once the pipeline is installed.

4) Modifying pipeline.py
------------------------

The actual workflow of the pipeline is defined in the pipeline.py. This script consists of three functions, initialize_all, 
run_all and report_all. The first one runs all initialize functions of all modules in the beginning, the second one controls the actual 
workflow and the final one calls all individual report functions in the end of a pipeline runs.

To begin with we need to import the star module at the top of the script:

.. code:: bash

	import rnaseq_pipeline.star

Then add star parameter initialization into the initialize_all:

.. code:: bash

	rnaseq_pipeline.star.init(param)


And finally we need to add the star module call. There is a submission function that makes our lives easier:

.. code:: bash

	HELPER.submit_job(param, pyfile, input_files, output_files, cores, mem_free)

Where:

``param``	Parameter object

``py_file``	specifies the entry point, as declared in setup.py ('``run_star``')

``input_files``	specifies the key in the parameter object, where the current working files are stored. Has to match up with the output_files key of the module that was run previously. 
	
``output_files``	Specifies the key to the parameter object in which to store the resulting files.

``cores``	(optional) Number of cores to be used.

``mem_free``	(optional) Free memory required on the node

.. code:: bash

      #do alignment if it's not just a fastqc run
        if not param['QC_and_trim_only']:
            if param['aligner'] == 'tophat':
                #running the aligner
                HELPER.submit_job(param,
                                  'run_tophat',
                                  input_files='fastq_files',
                                  output_files='bam_files',
                                  cores=param['qsub_num_processors'])
            if param['aligner'] == 'star':
                #running the aligner
                HELPER.submit_job(param,
                                  'run_star',
                                  input_files='fastq_files',
                                  output_files='bam_files',
                                  cores=param['qsub_num_processors'],
                                  mem_free='32G'))                    
            else:
                HELPER.writeLog('The selected aligner does not exist.', param)
                sys.exit(0)


In addition to the submit job function you can also write into the main log file using the HELPER.writeLog function.

5) Installing the changes
-------------------------

Once all the changes are done you can simply install and test them. For that you need the developer tools, 
if they are not already installed you can go into your github repository and install it using conda:

.. code:: bash

	cd CBMgithub/tools/RNASeq_pipeline
	conda install --file dev_requirements.txt
	
And once these are installed you can simply install the pipeline:

.. code:: bash

	python setup.py install

We already prepared unit test cases with only 10,000 reads that should allow you to test the pipeline in a reasonable amount
of time (depending on the cluster load usually <1h)

How to add a tool into the anaconda environment
-----------------------------------------------

A "star" package can be created by following the steps described in :ref:`rst_tutorial`


Supplement:
----------

**Building the Star index:**

A detailed description is provided under: 
`<https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf>`_.

Here is just a simple example on how we used to build our star index based on the human hg19 annotation. 

The manual suggests using gene annotations in .gtf format if they are available since it drastically increases 
its ability to find splice junctions correctly.

.. code:: bash

	STAR --runMode genomeGenerate \
		--genomeDir ~/annotation/star_index  \
		--genomeFastaFiles ~/annoation/genomes/genome.fa \
		--sjdbGTFfile ~/annotation/genes.gtf \
		--sjdbOverhang 100 \
		--runThreadN 1













