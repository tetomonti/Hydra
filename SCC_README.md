

## Setting Up Dependencies

The following workflow eliminates the direct use of most
modules. Basically you'll be setting up your own custom installs of each
dependency. To use this setup it's advised to not use the modules that
you're replacing, especially Python. To see what modules you have loaded,
use the `list` command:

    module list

To remove a package, use the `unload` command:

    module unload python2.7

Before we proceed, you'll need the `linga-proxy` module loaded if you're
working on `scc4.bu.edu`:

    module load linga-proxy

The python version you'll need to use is called `anaconda`:

    module load anaconda/2.1.0

This provides access to the `conda` tool, which is used to install
packages where ever you want. To setup your own install environment, use
the following command:

    conda create -p ./demo_env --copy python==2.7.3 pip --yes

This will create an install in the directory `./demo_env`, which means
it's setup relative to where you are currently working. It will install
python (currently pinned to the version used by the pipeline, but that can
change in the future) and pip, which you can use to install packages that
conda is not aware of. Basically try `conda install whatever` and if no
package is available, then try `pip install whatever`.

Once you've created the environment, then you can "activate" it, which
means make it the current one to use (this process is complementary to the
module system):

    source activate ./demo_env

Now install all the other dependencies:

    conda install --copy --yes \
       -c file:///restricted/projectnb/montilab-p/conda_channel \
       samtools bowtie2 cufflinks cutadapt fastqc\
       htseq pysam r==3.0.0 subread tophat

Traditionally conda is ment for hard to install python packages
(e.g. numpy, scipy and matplotlib), and conda is aware of an online
repository (hense the 'linga-proxy' requirement above). The montilab-p
project has a "channel" containing all the packages that the pipeline
uses, some that use R, some that use Python, all ready to be installed in
the same manor as a python package. The above command installs these
programs for you in your `./demo_env` directory.

Currently the pipeline is not configured to be installed in the same
manor, but that will change in the near future. When that happens, you can
install it along with all the other packages (in fact, you'll need to just
install the pipeline, and conda will automagically install all the
dependencies).


## Using conda-based setup

If you haven't already done this, activate your conda environment:

    module load anaconda
    source activate path/to/your/demo_env/can/be/relative

Now just run the pipeline:

    python path/to/CBMgithub/tools/RNASeq_pipeline/RNASeq_pipeline.py -p param.txt

If you ever want to stop using a conda environment, then 'deactivate' it:

    source deactivate

This will bring back the default anaconda install. 


## Using Git Repo

To get latest version:
git fetch origin

To add a file:
git add file (or *)

To push all the changes to the github repository:
git commit (filename or -a)
git push origin master

#nice and simple slide set that walks you through git + github
http://rogerdudler.github.io/git-guide/

Also there is a graphical interface: gitk