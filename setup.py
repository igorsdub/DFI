#! /usr/bin/python
"""Setuptools-based setup script for DFI.
For a basic installation just type the command::
  python setup.py install
"""

from setuptools import setup

setup(name='dfi',
      version='0.1.0',
      description='Dynamic Flexibility Index',
      author='Avishek Kumar',
      author_email='avishek.kumar@asu.edu',
      classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: POSIX',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Software Development :: Libraries :: Python Modules',
        ],
      packages=['dfi'],
      scripts=['./dfi/dfi_calc.py'],
      license='BSD',
      long_description=open('README.md').read(),
      install_requires=['numpy', 'pandas', 'pytest', 'matplotlib', 'biopython',
                        'scipy','seaborn']
      )
