.. _rst_tutorial:

=====================
 Conda Package Guide
=====================


Conda packages are created with `conda build`, which has a basic tutorial
`available online
<http://conda.pydata.org/docs/build_tutorials/pkgs.html>`_ .


Making a "star" Package
=======================

.. note::

   On SCC, you want to set the following environment variable to have your
   package automatically stored in a common location.

       export CONDA_BLD_PATH=`/restricted/projectnb/montilab-p/conda_packages/conda_build_space`

   Also, to access conda on scc, load the anaconda module

       module load anaconda/2.2.0


A conda package is built using two files. **build.sh** are the instructions you
would type manually to compile the package. **meta.yaml** is like a
configuration file, and it lists all the extra details needed to install
the packages, such as other requirements.

meta.yaml
---------

Copy one of the **meta.yaml** files for one of the packages found here::

  /restricted/projectnb/montilab-p/conda_packages

The first step is to set the package name and version in the **package**
section::

  package:
    name: star       
    version: "2.4.0j"
                   
Next you want to look for the download url for the source code. Usually
there's a link on the developer's site, or if it's a package built on SCC,
you can find a notes file in
/share/pkg/*package-name*/*version-number*/notes. This url and the file
name to be downloaded is needed in the **source** section::

  source:
    fn: STAR_2.4.0j.tar.gz
    url: https://github.com/alexdobin/STAR/archive/STAR_2.4.0j.tar.gz

Next you need to set some details in the **build** section, specifically
the *number* and *string* fields, which tags the build with a unique
name::

  build:
    number: 0       
    string: monti_0

Next you can use the **requirements** section to list requiremest for the
*build* process (sometimes you need special tools to build the program but
not run it) and the requirements for running the application. These fields
have lists as values, one item per line, preceeded by a dash. For example,
you might have::

  requirements:
    build:
      - cython
      - numpy
      - python

    run:
      - python
      - numpy

`star` has no requirements, so we don't need this section.

Next you can have a **test** section which you can use to define tests to
validate the build process.

Finally there is an **about** section that provides a link to the projects
*home* webpage and a description of the License agreement. 


build.sh
--------

This file is a shell script that lists all the steps you need to take to
build `star`. Look to the developer's website for help on the steps you
need to take for your application. When conda runs your **build.sh**
script it also sets up a bunch of environment variables to help the
process. The most important one is `$PREFIX`, which represents the
directory you should install your application to. Here's what the `star`
**build.sh** script looks like::

  cd source
  make
  mkdir $PREFIX/bin
  cp STAR $PREFIX/bin/
  ln -s $PREFIX/bin/STAR $PREFIX/bin/star

  
