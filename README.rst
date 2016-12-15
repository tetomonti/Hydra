
==================
 Hydra - a flexible RNASeq pipeline.
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

Hydra requires Java 1.6+. All other tools are wrapped in a Conda environment, see below. 


*1. Get Conda*


`conda` is available as part of the `miniconda package <http://conda.pydata.org/miniconda.html>`_.


.. note::

   If you're using a module system such as on the shared computing cluster (SCC) at Boston University you can just load a preinstalled version:

   .. code:: bash

     module purge
     module load anaconda2/4.2.0

*2. Set Environment and Build Paths*::

    export CONDA_ENVS_PATH=<WORKING_DIR>/conda-envs
    export CONDA_BLD_PATH=<WORKING_DIR>/conda-bld

*3. Create Dev Environment*

Use `conda` to install a basic developement environment::
  
  conda env create \
    reeder/dev_env \
    -p ./dev_env 

This will create a `./dev_env folder that stores all the files needed to run the pipeline. For more details, please refer to the Developer's Guide. (LINK)

*4. Install all necessary R packages*
Activate the environment so you have access to the R version that the pipeline uses::
 
  source activate ./dev_env

Start R::

  R

Within R install all necessary packages::

    require(devtools)
    install_git("http://github.com/nachocab/clickme.git")

Running the pipeline
====================

Activate the environment (follow the instructions provided by conda following the environment's
creation)::
 
  source activate ./dev_env
  
Then you can run the pipeline, while providing your parameter file.

Example::

   hydra-rnaseq -p param.txt


A detailed description of the parameter file is provided `here <docs/user_docs/parameters.rst>`__ and an example is located here `here <parameters_example.txt>`__ .


Additional documentation
========================

A detailed user documentation is provided `here <docs/user_docs/index.rst>`__. For advanced users that are contributing to the development such as adding additional modules a documentation is provided `here <docs/dev_docs/index.rst>`__


