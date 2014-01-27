#!/usr/bin/env python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

with open('README.md') as file:
    long_description = file.read()

setup(name="mip_family_analysis",
	version="0.5.2",
	author="Mans Magnusson",
	author_email="mans.magnusson@scilifelab.se",
	description=("A new tool for doing inheritance analysis and scoring in the mip pipeline."),
	long_description = long_description,
    packages={'Mip_Family_Analysis', 'Mip_Family_Analysis.Utils', 'Mip_Family_Analysis.Variants', 'Mip_Family_Analysis.Family', 'Mip_Family_Analysis.Models'},
    url='https://github.com/moonso/Mip_Family_Analysis',
    scripts=['scripts/run_mip_family_analysis.py'],
)
