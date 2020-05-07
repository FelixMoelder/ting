import sys
from setuptools import setup
# set __version__ and DESCRIPTION


with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='bio-ting',
    version='1.0.1',
    author='Felix Mölder',
    author_email='felix.moelder@uni-due.de',
    description='ting - T cell receptor interaction grouping',
    long_description=long_description,
    long_description_content_type="text/markdown",
    zip_safe=False,
    license='MIT License',
    url='https://github.com/FelixMoelder/ting',
    packages=['ting', 'scripts'],
    install_requires=['networkx>=2.4,<2.5', 'numpy>=1.17,<1.18', 'scipy>=1.3,<1.4'],
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
