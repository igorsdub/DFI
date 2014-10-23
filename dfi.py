#!/usr/bin/env python 
"""
USAGE: dfi.py PDBID

PDBID:    Protein DataBank Code 

"""

import sys 
import pdbmunge 
import pdbio 
import numpy as np 


if len(sys.argv) < 2:
    print __doc__ 
    exit()

if __name__ == "__main__":
    Verbose = True 
    gamma = 100
    

    #parse the input 
    pdbid = sys.argv[1] 
    
    #get the pdb, then extract the chain and then the alpha carbons and return the name of the file 
    fname = pdbmunge.extractA_CA(pdbid)
    ATOMS = [] 
    pdbio.pdb_reader(fname,ATOMS)

    x = []
    y = []
    z = [] 

    for atom in ATOMS:
        if(Verbose):
            print atom.atom_name
            print "%s %f %f %f"%(atom.atom_name,atom.x,atom.y,atom.z)
        if atom.atom_name == 'CA ':
            x.append(atom.x)
            y.append(atom.y)
            z.append(atom.z)
    
    x = np.array(x,dtype=float)
    y = np.array(y,dtype=float)
    z = np.array(z,dtype=float)

    #start computing the Hessian 
    numres = len(ATOMS)
    numresthree = 3 * numres
    
    hess = np.zeros((numresthree,numresthree))

    print hess 
    
