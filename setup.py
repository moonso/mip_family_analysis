#!/usr/bin/env python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


setup(name="Mip_Family_Analysis",
	version="0.3",
	author="Mans Magnusson",
	author_email="mans.magnusson@scilifelab.se",
	description=("A new tool for doing inheritance analysis and scoring in the mip pipeline."),
	long_description = open('README.md', 'r').read(),
    packages=['Mip_Family_Analysis', 'Mip_Family_Analysis.Utils', 'Mip_Family_Analysis.Variants', 'Mip_Family_Analysis.Family', 'Mip_Family_Analysis.Models'],
    scripts=['scripts/mip_family_analysis.py'],
)
