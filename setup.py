#! /usr/bin/python
"""Setuptools-based setup script for DFI.
For a basic installation just type the command::
  python setup.py install
"""

from setuptools import setup


def parse_requirements(fname):
    """
    Parses requirments.txt file into a list

    Parameters
    ----------
    fname: str
        name of the file with dependencies
        (e.g., requirements.txt)

    Returns
    -------
    dependencies: ls[str]
       ls of dependencies
    """
    with open(fname, 'r') as infile:
        dependencies = (
            infile.read().splitlines()
        )

    return dependencies


setup(name='dfi',
      version='0.10.0',
      description='Dynamic Flexibility Index',
      author='Avishek Kumar',
      author_email='avishek@asu.edu',
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
      scripts=['./dfi/dfi_calc.py',
               './dfi/uniprot_dfi.py',
               './dfi/fastaseq.py'],
      license='BSD',
      long_description=open('README.md').read(),
      install_requires=parse_requirements('requirements.txt')
      )
