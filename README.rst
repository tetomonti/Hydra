
==================
 RNASeq Pipeline.
==================

Features
========

* html reporting
* automatic parallelization
* automatic resume function
* extensive quality control
* interactive principal component analysis


Running the pipeline
====================


**Pre-requisites**

Below are the tools you need to install before running the pipeline. 

*1. Get Conda*


`conda` is available as part of the `miniconda package <http://conda.pydata.org/miniconda.html>`_.


.. note::

   If you're using a module system and conda is available as part of the
   `anaconda` module, you can install it with:

   .. code:: bash

     module load anaconda/2.2.0


*2. Create Dev Environment*

Use `conda` to install a basic developement environment::
  
  conda create \
    -p ./dev_env \
    -c 'file:///restricted/projectnb/montilab-p/conda_channel' \
    --yes \
    rnaseq_pipeline

This will create a `./dev_env folder that stores all the files needed to run the pipeline. For more details, please refer to the Developer's Guide.


**Running**

Activate the environment: follow the instructions provided by conda following the environment's
creation::
 
  source activate ./dev_env
  
Once you are done with your changes and installed them, you can run the pipeline. 

Example:

   RNASeq_pipeline_prototype.py -p param.txt
