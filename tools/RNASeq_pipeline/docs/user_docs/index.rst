
==============
 User's Guide
==============

"Pipeline name" is a semi-automated pipeline that runs a suite of tools required to perform analysis using high-throughput RNA Sequencing data. 

The pipeline takes list of user-specified parameters in the form of a parameter file (which also contains pointers to the pipeline scripts folder 
and raw files to be run using the pipeline) and does a complete run, with added flexibility in how programs are utilized 
(skipping certain steps, etc). 

An html report is generated, summarizing the output/results from each step, with links to the associated data (log files, QC reports etc.)

Table of Contents:

.. toctree::
   :maxdepth: 2
   
   parameters.rst
   running.rst
   output.rst


