#!/usr/bin/env python
"""
Download PDB
============

Description
-----------
Download pdb is suite of tools for interfacing
with the protein data bank (PDB)
"""
from __future__ import print_function
import numpy as np






def fetch_pdb(id, writetofile=True, Verbose=False):
    """
    Download pdb file and write out to file id.pdb

    Input
    ------
    id: str
       4 letter pdb code
    writetofile: bool
       Write out to file otherwise return file.
    Output
    ------
    id.pdb: file
       filename id.pdb
    """
    try:
        from StringIO import StringIO as stream
    except ImportError:
        from io import StringIO as stream
        from io import BytesIO as stream
    try:
        from urllib.request import urlopen
    except ImportError:
        from urllib2 import urlopen

    url = 'http://www.rcsb.org/pdb/files/%s.pdb' % id
    if(writetofile):
        with open(id + '.pdb', 'wb') as outfile:
            outfile.write(urlopen(url).read())
        if(Verbose):
            print("Wrote out %s.pdb" % (id))
    else:
        return stream(urlopen(url).read())
