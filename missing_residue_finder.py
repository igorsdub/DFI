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



def get_residues(pdbid,CAonly=True,noalc=True,chainA=False,chain_name='A',Verbose=False):
    "Returns the list of residues input is PDBCODE"
    
    url = 'http://www.rcsb.org/pdb/files/%s.pdb' %pdbid
    pdb = urllib.urlopen(url).readlines()
    residues = []
    for line in pdb:
            if line.startswith('ENDMDL'):
                print "MULTIPLE MODELS...USING MODEL1"
                return 
            
            if line.startswith('ATOM'):
                record = line[:6]
                atom_index = line[7:11]
              
                atom_name = line[13:16]
                if CAonly and not(atom_name=='CA '):
                    continue 
                
                alc = line[16] #alternate location
                if noalc and not((alc==' ' or alc=='A')):
                    continue 

                res_name = line[17:20]
                if(Verbose):
                    print res_name
                
                chainID=line[21]
                if chainA and not(chainID==chain_name):
                    continue 
                   
                res_index = line[22:27]
                if(Verbose):
                    print res_index 
                insert_code = line[26]
                x = line[31:38]
                y = line[39:46]
                z = line[47:54]
                occupancy = line[55:60]
                temp_factor = line[61:66]

                
                pdbnum = [pdbid, atom_name.strip(), res_name, chainID.strip(), int(res_index)]
                residues.append(pdbnum)
    if(Verbose):
        print residues[:10] 
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

def test_check_gaps_withgaps():
    assert check_gaps('1hcf') == '1hcf,TRUE,A', check_gaps('1hcf')

def test_check_gaps_input():
    assert check_gaps('1a4l') == '1a4l,FALSE,NA', check_gaps('1a4l')

if __name__ == "__main__":
    #input check 
    if len(sys.argv) < 2:
        print __doc__
        sys.exit()

    listfile = sys.argv[1]
    print "Input File: ", listfile 
    proc_list(listfile)
    












