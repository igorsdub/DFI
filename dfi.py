#!/usr/bin/env python 
"""
USAGE: dfi.py PDBID

PDBID:    Protein DataBank Code 

"""

import sys 

if len(sys.argv) < 2:
    print __doc__ 
    exit()

if __name__ == "__main__":
    Verbose = True 

    #parse the input 
    #get the pdb, then extract the chain and then the alpha carbons and return the name of the file 
    
    #read in the pdbfile 
