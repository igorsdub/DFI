#!/usr/bin/env python 
# coding: utf-8
"""
GAPFINDER
=========
 
Description
-------------
Find gaps and modify the fafsa sequence accordingly. 
- If the pdbstructure index jumps to a greater amount then increase the fafsa by the same amound 
- If the pdbstructure index jumps to a lower amount then renumber that part of the fafsa sequence. 

Usage
-------
```
gapfinder.py CSVFILE 
```
"""

import sys 
if len(sys.argv) < 2:
    print "No CSVFILE"
    print "------------------------"
    print __doc__
    sys.exit(1)

import pandas as pd
import fafsaseq 


def findgaps(gapfname):
    gapdata = pd.read_csv(gapfname,index_col='ResI')
    pdbind=gapdata.index.values 
    iterval=pdbind[0]
    
    for i in range(len(pdbind)):
    #print "pdbind",pdbind[i]
        if pdbind[i] == iterval:
            iterval += 1
        else:
            #print "GAP"
            if pdbind[i] > iterval:
                gapdist = pdbind[i] - iterval
            #print "gapdist",gapdist 
                gapdata.fafsa_ind.values[i:] += gapdist 
                iterval = gapdata.fafsa_ind.values[i] + 1
            #print "iterval before continue",iterval 
                continue
            elif pdbind[i] < iterval:
            #print "CASE B"
                gapdist = iterval - pdbind[i]
                gapdata.fafsa_ind.values[i:] -= gapdist 
                iterval = gapdata.fafsa_ind.values[i] + 1
           
    gapfnamels=gapfname.split('-')
    print gapfnamels
    outfilename=gapfnamels[0]+'-'+gapfnamels[1]+'-'+'nogaps'+'-'+gapfnamels[2]
    print outfilename
    gapdata.to_csv(outfilename)


if __name__ == "__main__":
    fname=sys.argv[1]

    #find the uniprot ids that go with this pdbid
    pdbid=fname.split('-')[0]
    print "PDBID: %s"%pdbid

    uniprols=fafsaseq.getuniprols(pdbid)
    print "Number of UNIPRO IDs: %d"%(len(uniprols))
    if len(uniprols) == 0:
        print "No UNIPROLS"
        exit() 

    gapfname_ls=fafsaseq.parsefafsaseq(fname,uniprols=uniprols)
    print gapfname_ls
    for gapfname in gapfname_ls:
        findgaps(gapfname)
