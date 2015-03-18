#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


readme = open('README.rst').read()
history = open('HISTORY.rst').read().replace('.. :changelog:', '')

requirements = [
    # TODO: put package requirements here
]

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='rnaseq_pipeline',
    version='1.0.0',
    description='A pipeline manager for RNASeq analysis',
    long_description=readme + '\n\n' + history,
    author='Monti Lab',
    author_email='smonti@bu.edu',
    url='https://github.com/montilab/CBMgithub',
    packages=[
        'rnaseq_pipeline',
    ],
    package_dir={'rnaseq_pipeline':
                 'rnaseq_pipeline'},
    include_package_data=True,
    install_requires=requirements,
    license="BSD",
    zip_safe=False,
    keywords='rnaseq_pipeline',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
    ],
    test_suite='tests',
    tests_require=test_requirements,
    scripts=['scripts/RNASeq_pipeline_prototype.py',
             'scripts/paired_ends_intersect.py',
             'scripts/run_bamqc.py'],
    entry_points={
        'console_scripts': [
            'bamqc.py=rnaseq_pipeline.bamqc:main',
            'cufflinks.py=rnaseq_pipeline.cufflinks:main',
            'cutadapt.py=rnaseq_pipeline.cutadapt:main',
            'fastqc.py=rnaseq_pipeline.fastqc:main',
            'featureCount.py=rnaseq_pipeline.featureCount:main',
            'htseq.py=rnaseq_pipeline.htseq:main',
            'matched_pairs.py=rnaseq_pipeline.matched_pairs:main',
            'tophat.py=rnaseq_pipeline.tophat:main',
            'trimmer.py=rnaseq_pipeline.trimmer:main']
    }
)
