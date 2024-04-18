#!/usr/bin/env python3

import subprocess
import sys
from setuptools import setup
from setuptools.command.build_ext import build_ext

class Build(build_ext):

    def run(self):
        protoc_command = ['make', 'all']
        if subprocess.call(protoc_command) != 0:
            sys.exit(-1)
        super().run()


setup(
    name='not-alike3',
    version='3.1.0',
    author='Javier Montalvo',
    author_email='buitrejma@gmail.com',
    py_modules=['not_alike3.nal', 'not_alike3.utils'],
    packages=['not_alike3'],
    python_requires='>=3.10',
    description='Not-alike3 commands pipeline finds unique dissimilar regions in a genome of interest by means of whole-genome systematic sequence comparison.',
    long_description = open('README.md', 'r').read(),
    install_requires =['click', 'pandas', 'biopython', 'numpy'],
    license = 'GNU General Public License v3 or later (GPLv3+)',
    url='https://www.github.com/exseivier/not-alike3-3.0.0',
#    cmdclass = {
#            'build_ext': Build
#        },
#    data_files = [
#            ('lib/python3.11/site-packages/not_alike3/biostruct', ['not_alike3/biostruct/libdnah.so'])
#        ],
#    include_package_data=True,
    entry_points = {
        'console_scripts' : [
                'not-alike3=not_alike3.nal:main'
                ]
            },
    classifiers = [
        'Development Status :: 2 - Alpha',
        'Environment :: Console',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.10',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
        )
