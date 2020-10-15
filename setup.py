'''
liv-utils (c) University of Liverpool. 2020

All rights reserved.

@author:  neilswainston
'''
# pylint: disable=exec-used
import setuptools


with open('README.md', 'r') as fle:
    _LONG_DESCRIPTION = fle.read()

with open('requirements.txt') as fle:
    _REQUIREMENTS = fle.read().splitlines()

setuptools.setup(
    name='liv_utils-utils',
    version='1.0',
    author='Neil Swainston',
    author_email='neil.swainston@liverpool.ac.uk',
    description='liv-utils',
    long_description=_LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    url='https://github.com/neilswainston/liv-utils',
    packages=setuptools.find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.7',
    install_requires=_REQUIREMENTS
)
