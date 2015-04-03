===============
 Documentation
===============

Introduction
============

To build the docs, from the root dirctory `CBMgithub/tools/RNASeq_pipeline`::

  make docs

To view the docs::

  firefox ./docs/_build/html/index.html

The docs are stored in `./docs/_build/html`, so copy that directory to an
HTTP server to host publicly. 

The directory `./docs/user_docs` contains documentation for using the
pipeline, where as the directory `./docs/dev_docs` contains the
documentation for developing code with the pipeline.

All documents are written using
`ReStructuredText <http://sphinx-doc.org/rest.html>`_ (.rst suffix). To
get the documents to show up in the table of contents, add them to the
corresponding `index.rst` files, e.g. `./docs/user_docs/index.rst`
contains a list of files for the user guide. 


Using Docstrings
================

When you use the make target `docs` from the root directory's Makefile,
you automatically run sphinx's api documentation script, and this creates
the appropriate \*.rst files in the `docs` directory. This automatically
makes a list of all modules and functions you write in the pipeline. The
docstrings are used as part of the documentation (`for example <https://pythonhosted.org/an_example_pypi_project/sphinx.html#function-definitions>`_). 

Basially, the unassigned string following a function's signature will be
used as part of the function's documentation::

  def add(x, y):
      """Return the sum of x and y"""
      return x + y

This will get converted into:

.. function:: add(x, y)

   Return the sum of x and y

You can add extra details to the docstring to encourage better
documentation formatting::

    def add(x, y):
      """Return the sum of x and y

      :param: x
      :param: y, must be same type as x
      :rtype: will be the same type as x and y
      """
      return x + y

.. function:: add(x, y)

   Return the sum of x and y
   
   :param: x
   :param: y, must be same type as x
   :rtype: will be the same type as x and y


Using sphinx markup (e.g. `:param:`) results in prettier html output, but
when reading the code, the result is quite ugly. The example linked to
above suggests using a format that wont
`automatically be formatted <https://pythonhosted.org/an_example_pypi_project/sphinx.html#full-code-example>`_ in sphinx, but is none-the-less organized. 
