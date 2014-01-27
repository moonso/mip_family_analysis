#!/usr/bin/env python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


setup(name="mip_family_analysis",
	version="0.4.1",
	author="Mans Magnusson",
	author_email="mans.magnusson@scilifelab.se",
	description=("A new tool for doing inheritance analysis and scoring in the mip pipeline."),
	long_description = open('README.md', 'r').read(),
    packages={'mip_family_analysis', 'mip_family_analysis.utils', 'mip_family_analysis.variants', 'mip_family_analysis.family', 'mip_family_analysis.models'},
    package_dir={'mip_family_analysis':'mip_family_analysis'},
    scripts=['scripts/run_mip_family_analysis.py'],
)
