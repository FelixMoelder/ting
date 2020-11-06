__author__ = "Felix Mölder"
__copyright__ = "Copyright 2020, Felix Mölder"
__email__ = "felix.moelder@uni-due.de"
__license__ = "MIT"

import sys

# set __version__ and DESCRIPTION

if sys.version_info < (3, 7):
    print("At least Python 3.7 is required.\n", file=sys.stderr)
    exit(1)

try:
    from setuptools import setup
except ImportError:
    print("Please install setuptools before installing ting.", file=sys.stderr)
    exit(1)

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='bio-ting',
    version='1.1.0',
    author='Felix Mölder',
    author_email='felix.moelder@uni-due.de',
    description='ting - T cell receptor interaction grouping',
    long_description=long_description,
    long_description_content_type="text/markdown",
    zip_safe=False,
    license='MIT License',
    url='https://github.com/FelixMoelder/ting',
    packages=['ting', 'scripts'],
    install_requires=['numpy>=1.17,<=1.19', 'scipy>=1.3,<=1.5'],
    entry_points={
        "console_scripts":
            ["ting=ting.ting:main", "imseq2ting=scripts.imseq2ting:main"]
        },
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]
)
