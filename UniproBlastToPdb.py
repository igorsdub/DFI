#!/usr/bin/env python 
# coding: utf-8
"""
UniproBlastToPDB
=================

Description 
------------
Given a Uniprot code will run blastp to retreive fasta sequence and blast xml output 
and parse the output in the following way:

Uniprot | PDBid | chain | query_to | query_from | Iter Query Len | e-value | Query Coverage | Sequence Identity 

Usage
-----

```
UniproBlastToPdb.py UNIPROid 
```

Output 
-------
- UniproID.fasta
- Unipro_blast.xml 
- Unipro.csv 
"""

import sys
if len(sys.argv) < 2:
    print __doc__
    sys.exit(1)

import random 
import Bio
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW
from Bio import SeqIO

from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Seq import Seq 


emails = ['avishek.kumar@asu.edu','akumar67@asu.edu','avishek.kumar@outlook.com','avishekkumar87@gmail.com','brandon.mac.butler@gmail.com']


def UniBLAST(code,Verbose=True):
    """
    Input
    ------
    Uniprot Code 
    UniBLAST(code)
    e.g. UniBLAST('O00238')
    
    Description
    -----------
    Outputs a Fasta sequence and
    Runs blasp looking through pdb database with the Uniprot code 
    
    Output
    ------
    - UniproID.fasta     FASTA Sequence 
    - UniproID_blast.xml Blast output in XML Format 
    """
    
    Entrez.email = random.choice(emails)
    if(Verbose):
        print "Using email: %s"%(Entrez.email)
    with open(code + ".fasta", "w") as out_file:
        net_handle = Entrez.efetch(db="nucleotide", id=code, rettype="fasta")
        out_file.write(net_handle.read())
    
    print "Running blastp"
    result_handle = NCBIWWW.qblast("blastp", "pdb", code)
    print "Done running blastp"
    with open(code + "_blast.xml", "w") as save_file:
        save_file.write(result_handle.read())
    result_handle.close()

def getfasta(code,Verbose=True):
    """
    Input
    -----
    Uniprot code 

    Description
    ------------
    getfasta(code,Verbose=True)
    Query Webserver for FASTA sequence 
    
    Output
    ------
    code.fasta FASTA sequence 
    """
    from Bio import Entrez 
    Entrez.email = random.choice(emails)
    if(Verbose):
        print "Using email: %s"%(Entrez.email)
    with open(code + ".fasta", "w") as out_file:
        net_handle = Entrez.efetch(db="nucleotide", id=code, rettype="fasta")
        out_file.write(net_handle.read())
    

def parseBlastFile(xmlfil):
    """
    Input 
    -----
    Uniprot XML Output  
    parseBlastFile(xmlfil)
    e.g. parsebBlastFile('O00238_blast.xml')
    
    Description
    -----------
    Parses the following information out of the bast output xml file. 
    
    Uniprot | PDBid | chain | query_to | query_from | Iter Query Len | e-value | Query Coverage | Sequence Identity 
    
    Output 
    -------
    Uniprot.csv 
    """
    NP_id = xmlfil.split('_')[0]
    
    result_handle = open(xmlfil)
    blast_record = NCBIXML.read(result_handle)
    result_handle.close()
    
    #E_VALUE_THRESH = 1E-25
    
    outfilname = NP_id+'.csv'
    
    with open(outfilname,"w") as out_file:
        print "Writing output to %s"%(outfilname) 
        out_file.write('Uniprot,PDBid,chain,query_to,query_from,IterQueryLen,e-value,QueryCov,SeqId\n')
        sequencequeryLength = blast_record.query_length
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                first = float(hsp.identities)
                second = len(hsp.query)
                identity = 100*float(first/second)
                identity = int(round(identity,0))
                coverage =  round(100*float(hsp.query_end - hsp.query_start)/sequencequeryLength,0) 
                             
                line1=alignment.title
                b=line1.split('|')
                pdbid = str(b[3])
                
                out_file.write(NP_id+",")
                out_file.write(pdbid+",")
                line2=b[4]
                chain=line2.split()
                               
                out_file.write(str(chain[0]) +",")
                out_file.write(str(hsp.query_end)+",")
                out_file.write(str(hsp.query_start)+",")
                out_file.write(str(sequencequeryLength)+",")
                out_file.write(str(hsp.expect) +",")
                out_file.write("%f"%coverage +",")
                out_file.write(str(identity) +"\n")

    
if __name__ == "__main__":
    blastfile = sys.argv[1]
    #UniBLAST(code)  
    parseBlastFile(blastfile)

