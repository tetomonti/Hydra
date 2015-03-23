#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


readme = open('README.rst').read()
history = open('HISTORY.rst').read().replace('.. :changelog:', '')

with open('requirements.txt') as fin:
    requirements = map(line.strip() for line in fin.readlines())

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
            'run_bamqc=rnaseq_pipeline.bamqc:main',
            'run_cufflinks=rnaseq_pipeline.cufflinks:main',
            'run_cutadapt=rnaseq_pipeline.cutadapt:main',
            'run_fastqc=rnaseq_pipeline.fastqc:main',
            'run_featureCount=rnaseq_pipeline.featureCount:main',
            'run_htseq=rnaseq_pipeline.htseq:main',
            'run_matched_pairs=rnaseq_pipeline.matched_pairs:main',
            'run_tophat=rnaseq_pipeline.tophat:main',
            'run_trimmer=rnaseq_pipeline.trimmer:main']
    }
)
