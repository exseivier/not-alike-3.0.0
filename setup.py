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
    name='not-alike',
    version='3.0.0',
    author='Javier Montalvo',
    author_email='buitrejma@gmail.com',
    py_modules=['not_alike.nal', 'not_alike.utils'],
    packages=['not_alike'],
    python_requires='>=3.9',
    description='Command pipeline that finds not alike regions of query genome compared to at least one genome form a huge list of different genomes.',
    long_description = open('README.md', 'r').read(),
    install_requires =['click', 'pandas', 'biopython'],
    license = 'GNU General Public License v3 or later (GPLv3+)',
    url='https://www.github.com/exseivier/not-alike',
    cmdclass = {
            'build_ext': Build
        },
    data_files = [
            ('lib/python3.11/site-packages/not_alike/biostruct', ['not_alike/biostruct/libdnah.so'])
        ],
    include_package_data=True,
    entry_points = {
        'console_scripts' : [
                'not-alike=not_alike.nal:main'
                ]
            },
    classifiers = [
        'Development Status :: 1 - Pre-Alpha',
        'Environment :: Console',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
        )
