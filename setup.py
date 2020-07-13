'''
genegeniebio-utils (c) GeneGenie Bioinformatics Ltd. 2020

All rights reserved.

@author:  neilswainston
'''
# pylint: disable=exec-used
import setuptools

with open('README.md', 'r') as fh:
    _LONG_DESCRIPTION = fh.read()

setuptools.setup(
    name='genegeniebio-utils',
    version='1.01',
    author='Neil Swainston',
    author_email='neil.swainston@genegenie.bio',
    description='genegeniebio-utils',
    long_description=_LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    url='https://github.com/genegeniebio/genegeniebio-utils',
    packages=setuptools.find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.7',
)
