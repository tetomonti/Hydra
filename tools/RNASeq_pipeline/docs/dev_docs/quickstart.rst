
=============
 Quick Start
=============

Below are all the steps you need to make changes to the pipeline. 

1. Get Conda
============

`conda` is available as part of the `miniconda package <http://conda.pydata.org/miniconda.html>`_.


.. note::

   If you're using SCC (BU's cluster), conda is available as part of the
   `anaconda` module

   .. code:: bash

     module load anaconda/2.2.0


2. Create Dev Environment
=========================

Use `conda` to install a basic developement environment::
  
  conda create \
    -p ./dev_env \
    -c 'file:///restricted/projectnb/montilab-p/conda_channel' \
    --yes \
    rnaseq_pipeline

Alternatively, you can just install the dependencies:

  conda create \
    -p ./dev_env \
    -c 'file:///restricted/projectnb/montilab-p/conda_channel' \
    --yes \
    samtools \
    bowtie2 \
    cufflinks \
    cutadapt \
    fastqc \
    htseq \
    pysam \
    r=3.0.0 \
    subread \
    tophat \
    matplotlib \
    numpy \
    python=2.7*

 And then you can install the package from source (see below).


What just happend?
------------------

1. `-p ./dev_env` -- use `./dev_env` as the installation path.

   This will create a `./dev_env/bin`, `./dev_env/lib`, etc directories to
   store all the files needed to run the programs you install using
   `conda`.

2. `-c 'file:///restricted/projectnb/montilab-p/conda_channel'` -- use a
   custion "channel".

   When `conda` downloads the files needed to install your programs, it
   looks in certain predefined locations. The `-c` flag adds a location to
   this search path, and it's the way to add custon packages to
   install. The Monti lab has a custom channel on SCC, which you've
   included in the search path with the above instruction.

3. `--yes` -- say "yes" to all posed questions.

   `conda` checks in with you before it actually does anything, and this
   instruction saves a bit of time. You can walk away from the computer,
   and when you get back, everything should just work.

4. `rnqseq_pipeline` -- the package in question.

   The channel configured above (with `-c` option) has a package called
   `rnaseq_pipeline` that contains the python and R code for the
   pipeline. This package was created with :doc:`conda build
   <making_conda_packages>`.

3. Activate the Environment
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

4. Get the Source
=================

Use git to clone the repository (using the group's repo or your own github
copy)::

  git clone https://github.com/montilab/CBMgithub.git


To get the current branch use::

  git clone https://user@github.com/montilab/CBMgithub.git -b v2.0.0



5. Install Developer Tools
==========================

The source contains a "dev_requirements.txt" file that lists all the
packages used in development. Install these using conda::

  cd CBMgithub/tools/RNASeq_pipeline
  conda install --file dev_requirements.txt

.. warning::

   Make sure you're installing into the correct environment. `which
   python` will print the path to the python being used, and make sure
   it's the one in the development environment. If not, see
   `3. Activate the Environement`_ to setup the environment. 


6. Make Your Edits and Install
==============================

Once you've made your edits, you can install them as you would any
standard python package::

  python setup.py install


7. Test whether your changes broke the pipeline
===============================================

Once you are done with your changes and installed them, try running one 
of the toy examples to make sure there were no unintended side effects::

   cd /restricted/projectnb/montilab-p/projects/pipeline_dev/unit_tests/human_paired_end
   RNASeq_pipeline_prototype.py -p param.txt


Once you're sure everything works, use git to commit them::

   git add <whatever_file_you_changed>
   git commit -m "Add a meaningful commit message"
   git push


To see the files that have been changed and need to be committed use::
   git status

