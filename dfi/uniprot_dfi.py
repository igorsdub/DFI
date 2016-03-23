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
import sys 
import UniproBlastToPdb as uni
from download_pdb import fetch_pdb 
import dfi_calc 

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
        blastfile = code+'_blast.xml'
        csvfile = code+'.csv'
        print "Blasting"
        uni.UniBLAST(code)  
        print "ParseFile"
        uni.parseBlastFile(blastfile)
        print "Get top hit"
        pdbid = uni._gettophit(csvfile)
        if pdbid == None:
            print code, "No Good PDB hit"
        else:
            print "Top Hit: ",pdbid 
        outfile = code+'_'+pdbid+'-dfianalysis.csv' 
        print "Running DFI"
        fetch_pdb(pdbid,Verbose=True)
        pdbfile = pdbid+'.pdb'
        dfi_calc.calc_dfi(pdbfile,pdbid,dfianalfile=outfile)
    

if __name__ == "__main__" and len(sys.argv) < 2:
    print __doc__
else:
    uniprotcodes = sys.argv[1:]
    uniproDFI(uniprotcodes)
    

    
