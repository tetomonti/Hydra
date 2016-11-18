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
    packages=['hydra_pkg'],
    include_package_data=True,
    install_requires=requirements,
    license="BSD",
    zip_safe=False,
    keywords='hydra',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache License',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
    ],
    test_suite='tests',
    tests_require=test_requirements,
    scripts=['scripts/hydra-rnaseq',
             'scripts/paired_ends_intersect.py',
             'scripts/run_bamqc.py'],
    entry_points={
        'console_scripts': [
            'run_bamqc=hydra_pkg.bamqc:main',
            'run_cufflinks=hydra_pkg.cufflinks:main',
            'run_cutadapt=hydra_pkg.cutadapt:main',
            'run_fastqc=hydra_pkg.fastqc:main',
            'run_featureCount=hydra_pkg.featureCount:main',
            'run_htseq=hydra_pkg.htseq:main',
            'run_matched_pairs=hydra_pkg.matched_pairs:main',
            'run_tophat=hydra_pkg.tophat:main',
            'run_star=hydra_pkg.star:main',
            'run_bowtie2=hydra_pkg.bowtie2:main']
    }
)
