"""
Download PDB
============

Description
-----------
Download pdb is suite of tools for interfacing
with the protein data bank (PDB)
"""
from __future__ import print_function


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
    from StringIO import StringIO
    import urllib
    url = 'http://www.rcsb.org/pdb/files/%s.pdb' % id
    if(writetofile):
        with open(id + '.pdb', 'w') as outfile:
            outfile.write(urllib.urlopen(url).read())
        if(Verbose):
            print("Wrote out %s.pdb" % (id))
    else:
        return StringIO(urllib.urlopen(url).read())
