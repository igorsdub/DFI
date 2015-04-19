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
    uniproid=''
   
    if uniprotls:
        print uniproURL+uniproid+'.fasta'
        response = urllib2.urlopen(uniproURL+uniproid+'.fasta')
    else:
        print "Taking from the PDB"
        response = urllib2.urlopen(pdbURL+pdbname)
  
    html = response.read()    
    html = html.split('\n')
    fseq=''
    for line in html:
        print line
        if line.startswith('>'):
            continue
        fseq = fseq + line

    numseq=4
    for i in range(len(fseq)):
        print fseq[i:i+numseq], smallseq[:numseq]
        if smallseq[:numseq] == fseq[i:i+numseq]:
            match = i+1
            print "Found a match at %d"%(i+1)
            break

    data['fafsa']=pd.Series( range(match,len(smallseq)+match), index=data.index)
    if uniprotls:
        outfile=pdbname+'-dfianalysis.csv'
    else:
        outfile=pdbname+'-'+uniproid+'-dfianalysis.csv'
    data.to_csv(outfile)
    print "Writing out to: " + outfile 

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print __doc__
        print sys.exit()
    parsefafsaseq(sys.argv[1],uniprot=False)
