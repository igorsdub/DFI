#!/usr/bin/env python 
"""
FAFSA SEQ
=========

Description
------------
Downloads the fafsa sequence from the pdb and sees where the sequence begins in
the fafsa sequence and adds it to the csv file. 

Usage
-----
```
fafsaseq.py DFICSVFILE
```
"""

import pandas as pd
import numpy as np 

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


def parsefafsaurl(html):
    "Parses html fafasa sequence and returns a string of the sequence"
    html = html.split('\n')
    fseq=''
    for line in html:
        print line
        if line.startswith('>'):
            continue
        fseq = fseq + line
    return fseq 

def compareseq(smallseq,fseq,numseq=4):
    "Compare sequence to find contiguous sequence"
    for i in range(len(fseq)):
        print fseq[i:i+numseq], smallseq[:numseq]
        if smallseq[:numseq] == fseq[i:i+numseq]:
            match = i+1
            print "Found a match at %d"%(i+1)
            return match 
    print "No match"
    return False 


def parsefafsaseq(fname,uniprols=None):
    """
    Parse the fafas seq using the csv filename and uniprotids
    """
    
    data = pd.read_csv(fname,index_col='ResI')
    
    #convert three letter code to one letter code sequence 
    bigseq=np.array(data['Res'].values,dtype=str)
    smallseq=''
    for threeseq in bigseq:
        smallseq = smallseq + str(mapres[threeseq])

    import urllib2
    pdbURL="http://www.pdb.org/pdb/files/fasta.txt?structureIdList="
    uniproURL="http://www.uniprot.org/uniprot/"
   
    pdbname=fname.split('-')[0]
       
    if uniprols:
        for uniproid in uniprols:
            print uniproid
            print uniproURL+uniproid+'.fasta'
            response = urllib2.urlopen(uniproURL+uniproid+'.fasta')
            fseq = parsefafsaurl(response.read())
            match = compareseq(smallseq,fseq,numseq=4)
            if match:
                ind_match=match - 1 
                matchseq = [f for f in fseq[ind_match:ind_match+len(smallseq)+ind_match] ]
                print fseq 
                print "---"
                print smallseq 
                print len(matchseq)
                print len(smallseq)
                if len(matchseq) < len(smallseq):
                    ntimes = len(smallseq) - len(matchseq) 
                    for i in range(ntimes):
                        matchseq.append('NA')
                
                print "fafalen:" + str(len(fseq))
                data['fafsa_seq']=pd.Series( matchseq, index=data.index)
                data['fafsa_ind']=pd.Series( range(match,len(smallseq)+match), index=data.index)
                data['unipro'] = pd.Series( [ uniproid for i in range(match, len(smallseq)+match) ], index=data.index)
                outfile=pdbname+'-'+uniproid+'-dfianalysis.csv'
                print "Writing out to: " + outfile
                data.to_csv(outfile)
    else:
        print "Taking from the PDB"
        response = urllib2.urlopen(pdbURL+pdbname)
        fseq = parsefafsaurl(response.read())
        match = compareseq(smallseq,fseq,numseq=4)
        if match: 
            data['fafsa_ind']=pd.Series( range(match,len(smallseq)+match), index=data.index)
            outfile=pdbname+'-dfianalysis.csv'
            data.to_csv(outfile)
            

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print __doc__
        print sys.exit()
    parsefafsaseq(sys.argv[1],uniprols=None)
