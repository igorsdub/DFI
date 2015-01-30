#!/usr/bin/env python 
"""
=========
FAFSA SEQ
=========

Description
------------
Downloads the fafsa sequence from the pdb and sees where the sequence begins in
the fafsa sequence and adds it to the csv file. 

Usage
-----
fafsaseq.py DFICSVFILE
"""

import sys

if len(sys.argv) < 2:
    print __doc__
    print sys.exit()

import pandas as pd
import numpy as np 

#mapres={}
#with open('RESdict.txt','r') as infile:
#    for line in infile:
#        ls_line=line.strip('\n').split(':')
#        mapres[ls_line[0]] = ls_line[1]

mapres={'ALA':'A',
'CYS':'C',
'ASP':'D',
'GLU':'E',
'PHE':'F',
'GLY':'G',
'HIS':'H',
'ILE':'I',
'LYS':'K',
'LEU':'L',
'MET':'M',
'PRO':'P',
'ARG':'R',
'GLN':'Q',
'ASN':'N',
'SER':'S',
'THR':'T',
'TRP':'W',
'TYR':'Y',
'VAL':'V'}

print mapres

#fname='1iau-dfianalysis.csv'
fname=sys.argv[1]
data = pd.read_csv(fname,index_col='ResI')

bigseq=np.array(data['Res'].values,dtype=str)
smallseq=''
for threeseq in bigseq:
    smallseq = smallseq + str(mapres[threeseq])

import urllib2
URL="http://www.pdb.org/pdb/files/fasta.txt?structureIdList="
#pdbname="1IAU"
pdbname=fname.split('-')[0]
response = urllib2.urlopen(URL+pdbname)
html = response.read()    
html = html.split('\n')
fseq=''
for line in html:
    print line
    if line.startswith('>'):
        continue
    fseq = fseq + line

for i in range(len(fseq)):
    print fseq[i:i+10], smallseq[:10]
    if smallseq[:10] == fseq[i:i+10]:
        match = i+1
        print "Found a match at %d"%(i+1)
        break

data['fafsa']=pd.Series( range(match,len(smallseq)+match), index=data.index)

data.to_csv(fname)
