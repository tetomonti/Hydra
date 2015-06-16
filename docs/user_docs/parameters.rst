===================================
Parameters and files specifications
===================================

Phenotype file
==============

The phenotype file is a separate tab delimited file with the following specifications:

- First row contains the header information
- The first column contains the location of the raw file (for paired end sequencing data the first 2 columns)
- The second column (third column for paired end) contains the sample names, which are used for output files
- The sample names should be called *sample_names* in the header, otherwise the pipeline won't automatically create an R/Bioconductor ExpressionSet
- In addition to raw file location and sample name the file can also contain any number of additional phenotype columns
- The last column is used for a PCA that is automatically created using the FPKM normalized data

.. note::

   All the file must be tab delimited, included names of the columns in the header row

Example::

	FastQ_files	sample_name	genotype	pregnenolone_treated	Group

	Sample_Ctr-MO_D1.R1.fastq.gz	Sample_Ctr-MO_D1.R1	control	no	Ctrl_DMSO

	Sample_Ctr-MO_D2.R1.fastq.gz	Sample_Ctr-MO_D2.R1	control	no	Ctrl_DMSO

	Sample_Ctr-MO_D3.R1.fastq.gz	Sample_Ctr-MO_D3.R1	control	no	Ctrl_DMSO

	Sample_Ctr-MO_D4.R1.fastq.gz	Sample_Ctr-MO_D4.R1	control	no	Ctrl_DMSO

	Sample_Ctr-MO_PN1.R1.fastq.gz	Sample_Ctr-MO_PN1.R1	control	yes	Ctrl_PN


Parameter file
==============

All other parameters are specified in the parameter file. This file follows the format ``<KEY> = <VALUE>``, with all test after a '#' removed. 
It does not allow spaces within keys or values. An example is available (parameters_example.txt) including a description for each parameter. 
The following list explains some of the more important parameters:

``working_dir``: specifies the location where the intermediate results, deliverables and final report will be output. 

``raw_filenames``: file with location of raw files and phenotype information as described above

``paired``: Specifies if the data are paired end

``clean_run``: Indicates if the pipeline should be run in resume mode or from scratch. If it runs in resume mode it will only run the steps of the pipeline that haven't finished

``verbose``: output the main log also onto the console

``aligner``: tophat (standard) or skip if one wants to start with bam files, which are then specified in the raw file directory. 

``QC_and_trim_only``: runs only the initial preprocessing step to make sure the adapter clipping was successful before continuing. 
One you confirm the QC looks alright you can run the pipeline in resume mode, which the skips the fastqc runs and adapter clipping.

	
There are also parameters for every single module (most important are the genome and annotation files) that have to be specified correctly. 

Annotation files can be downloaded at `<http://useast.ensembl.org/Homo_sapiens/Info/Index?redirect=no>`_ 


