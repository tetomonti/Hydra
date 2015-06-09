
=============
 Quick Start
=============

Below are all the steps you need to make changes to the pipeline. 

1. Get Conda
============

`conda` is available as part of the `miniconda package <http://conda.pydata.org/miniconda.html>`_.


.. note::

   If you're using a module system such as on the shared computing cluster (SCC) at Boston University you can just load a preinstalled version:

   .. code:: bash

     module purge
     module load anaconda/2.2.0


2. Create Dev Environment
=========================
Use `conda` to install a basic developement environment::
  
  conda create \
    -p ./dev_env \
    -c https://conda.binstar.org/montilab \
    --yes \
    hydra

This will create a `./dev_env folder that stores all the files needed to run the pipeline. For more details, please refer to the Developer's Guide. (LINK)

3. Install all necessary R packages
===================================
Activate the environment so you have access to the R version that the pipeline uses::
 
  source activate ./dev_env

Start R::

  R

Within R install all necessary packages::

  source("http://bioconductor.org/biocLite.R")
  biocLite(c("Biobase",'edgeR'),ask=F)
  install.packages(c("devtools","knitr", "yaml", "rjson"), repos='http://cran.us.r-project.org')
  devtools::install_github("clickme", "nachocab")

4. Activate the Environment
============================

Follow the instructions provided by conda following the environment's
creation::

  source activate ./dev_env

This essentially adds `./dev_env/bin` to your `$PATH` environment
variable, making whatever that's installed there available to run in the
current session, including python, bowtie2, etc. 

.. note::

   You need to do this each time (login session) you want to use this
   installation of the pipeline.

5. Get the Source
=================

Use git to clone the repository (using the group's repo or your own github
copy)::

  git clone https://github.com/montilab/Hydra.git


To get the current branch use::

  git clone https://user@github.com/montilab/Hydra.git -b v2.0.0



6. Install Developer Tools
==========================

The source contains a "dev_requirements.txt" file that lists all the
packages used in development. Install these using conda::

  conda install -r dev_requirements.txt

This will force you to downgrade the readline module however you can just install this single module again::

  conda install \
     --override-channels \
     -c https://conda.binstar.org/montilab \
     readline-6.2.3


.. warning::

   Make sure you're installing into the correct environment. `which
   python` will print the path to the python being used, and make sure
   it's the one in the development environment. If not, see
   `4. Activate the Environement`_ to setup the environment.


7. Make Your Edits and Install
==============================

Once you've made your edits, you can install them as you would any
standard python package, by switching to the repository::

  python setup.py install


8. Test whether your changes broke the pipeline
===============================================

Once you are done with your changes and installed them, try running one 
of the toy examples to make sure there were no unintended side effects::

   cd /restricted/projectnb/montilab-p/projects/pipeline_dev/unit_tests/human_paired_end
   hydra -p param.txt


Once you're sure everything works, use git to commit them::

   git add <whatever_file_you_changed>
   git commit -m "Add a meaningful commit message"
   git push


To see the files that have been changed and need to be committed use::

   git status

