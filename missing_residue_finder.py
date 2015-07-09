#!/usr/bin/env python 
# coding: utf-8
"""
Missing Residue Finder 
======================

Description
-----------
Checks if there are are missing residues by finding gaps in the numbering. 


Usage
-----
./missing_residue_finder.py LIST
"""

import numpy as np 
import pandas as pd
import urllib
import sys 



def get_residues(pdbid):
    "Returns the list of residues input is PDBCODE"
    
    url = 'http://www.rcsb.org/pdb/files/%s.pdb' %pdbid
    pdb = urllib.urlopen(url).readlines()
    residues = []
    for line in pdb:
        list = line.split()
        id = list[0]
        if id == 'ATOM':
            type = list[2]
            if type == 'CA':
                pdbnum = [pdbid, type, list[3], list[4], list[5]]
                residues.append(pdbnum)
    return residues




def get_chains(residues):
    "Returns a list of chains"
    chains = []
    for res in residues:
        chains.append(res[3])
    chains = np.unique(chains)
    return chains




def get_resnums_chain(residues,chainID,Verbose=False):
    "get residue numbers for a given chain"
    if(Verbose):
        print chainID
    resnums = []
    for line in residues:
        if (line[3] == chainID):
            resnums.append(line[4])
    return resnums


def check_gaps(pdbid):
    "check for gaps based on pdbID"
    residues = get_residues(pdbid)
    chains = get_chains(residues)
    for chain in chains:
        resnums = get_resnums_chain(residues,chain)
        startres = int(resnums[0])
        endres = int(resnums[-1])
        nogaps = [x for x in range(startres,endres+1)]
        if (len(resnums) != len(nogaps)):
            return "{pdb},{val},{chain}".format(pdb=pdbid,val='TRUE',chain=chain)
        
    return "{pdb},{val},{chain}".format(pdb=pdbid,val='FALSE',chain='NA')
    

def proc_list(listfil):
    with open(listfil,'r') as infile:
        for pdbid in infile:
            print check_gaps(pdbid.strip())

def test_check_gaps():
    assert check_gaps('1hcf') == '1hcf,TRUE,A'


if __name__ == "__main__":
    #input check 
    if len(sys.argv) < 2:
        print __doc__
        sys.exit()

    listfile = sys.argv[1]
    print "Input File: ", listfile 
    proc_list(listfile)
    












