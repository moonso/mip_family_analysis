#!/usr/bin/env python

from setuptools import setup, find_packages, Command


setup(name="Mip_Family_Analysis",
	version="0.1",
	author="Mans Magnusson",
	author_email="mans.magnusson@scilifelab.se",
	description=("A new tool for doing inheritance analysis and scoring in the mip pipeline."),
	long_description = open('README.md', 'r').read(),
    packages={'Mip_Family_Analysis', 'Mip_Family_Analysis.Utils', 'Mip_Family_Analysis.Variants', 'Mip_Family_Analysis.Family', 'Mip_Family_Analysis.Models'},
	package_dir={'Mip_Family_Analysis': 'Mip_Family_Analysis'},
    scripts=['scripts/mip_family_analysis.py'],
)
