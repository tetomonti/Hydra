#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


readme = open('README.rst').read()
history = open('HISTORY.rst').read().replace('.. :changelog:', '')

with open('requirements.txt') as fin:
    requirements = [line.strip() for line in fin.readlines()]

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='hydra',
    version='2.0.0',
    description='A pipeline manager for RNASeq analysis',
    long_description=readme + '\n\n' + history,
    author='Monti Lab',
    author_email='smonti@bu.edu',
    url='https://github.com/montilab/Hydra',
    packages=[
        'hydra',
    ],
    package_dir={'hydra':
                 'hydra'},
    include_package_data=True,
    install_requires=requirements,
    license="BSD",
    zip_safe=False,
    keywords='hydra',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
    ],
    test_suite='tests',
    tests_require=test_requirements,
    scripts=['scripts/hydra.py',
             'scripts/paired_ends_intersect.py',
             'scripts/run_bamqc.py'],
    entry_points={
        'console_scripts': [
            'run_bamqc=hydra.bamqc:main',
            'run_cufflinks=hydra.cufflinks:main',
            'run_cutadapt=hydra.cutadapt:main',
            'run_fastqc=hydra.fastqc:main',
            'run_featureCount=hydra.featureCount:main',
            'run_htseq=hydra.htseq:main',
            'run_matched_pairs=hydra.matched_pairs:main',
            'run_tophat=hydra.tophat:main',
            'run_star=hydra.star:main']
    }
)
