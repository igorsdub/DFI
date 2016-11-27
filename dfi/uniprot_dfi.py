#!/usr/bin/env python
"""
===
DFI
===

Description
-----------

Given a list of uniprot IDs:
dfi.py will find do a blast search on the NCBI
to find the highest hit PDB and calculate the DFI profile
of that pdb

Example
--------
./dfi.py P42771

"""
from __future__ import print_function
import sys
import dfi.UniproBlastToPdb as uni


def uniproDFI(uniprotcodes):
    """
    uniproDFI take a list of uniprot codes,
    finds the top pdb hit and then computes
    the dfi profile

    Input
    -----
    uniprotcodes: ls
       ls of uniprot codes to run DFI on

    """
    for code in uniprotcodes:
        blastfile = code + '_blast.xml'
        print("Blasting")
        uni.UniBLAST(code)
        print("ParseFile")
        uni.parseBlastFile(blastfile)


if __name__ == "__main__" and len(sys.argv) < 2:
    print(__doc__)
else:
    uniprotcodes = sys.argv[1:]
    uniproDFI(uniprotcodes)
