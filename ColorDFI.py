#!/usr/bin/env python 
# coding: utf-8
"""
Color by DFI
=============

Description
-------------
This program should read in a pdb and the dfi analysis of the pdb and then color
the CA accordingly to the dfi. 

Usage
------
```
CMDLINE:ColorDFI.py CSVFIL PDBFIL
colorbydfi(CSVFIL,PDBFIL,Verbose)

``` 

- CSVFIL: DFI CSV FILE 
- PDBFIL: Corresponding PDBFIL 
- Verboose: Boolean for debugging 
"""

import sys 


def colorbydfi(CSVFIL,PDBFIL,Verbose=False,colorbyparam='pctdfi',outfile=None):
    """
    Color by DFI
    =============

    Description
    -------------
    This program should read in a pdb and the dfi analysis of the pdb and then color
    the CA accordingly to the dfi. 

    Usage
    ------
    ```
    CMDLINE:ColorDFI.py CSVFIL PDBFIL
    colorbydfi(CSVFIL,PDBFIL)

    ``` 

    - CSVFIL: DFI CSV FILE 
    - PDBFIL: Corresponding PDBFIL 
    - Verbose: Boolean for debugging 
    """
    import pandas as pd
    import pdbio as io
    
    if type(CSVFIL) == str:
        data = pd.read_csv(CSVFIL)
        pdbid = CSVFIL.split('-')[0]
    else:
        data = CSVFIL
        pdbid = outfile.split('-')[0]
   
    if(Verbose):
        print "CSVFIL: %s"%(CSVFIL)
        print "PDBFIL: %s"%(PDBFIL)
        print "pdbid: %s"%(pdbid)
        print data[:10]
        print "Reading in: %s"%(CSVFIL)

    ATOMS = []
    io.pdb_reader(PDBFIL,ATOMS)
        
    for i in range(len(ATOMS)):
        if True:
            resind = ATOMS[i].res_index 
            chainind = ATOMS[i].chainID
            val= data[ ( data.ResI == int(resind) ) & ( data.ChainID == chainind ) ][colorbyparam].values[0]
            if Verbose:
                print ATOMS[i].res_index, ATOMS[i].temp_factor,val
            ATOMS[i].temp_factor = val
        else:
            ATOMS[i].temp_factor = 0. 
    if(outfile):
        io.pdb_writer(ATOMS,filename=outfile)
    else:
        io.pdb_writer(ATOMS,filename=pdbid+'-dficolor.pdb')

if __name__ == "__main__" and len(sys.argv) < 2:
    print __doc__
    exit()

if __name__ == "__main__":
    colorbydfi(sys.argv[1],sys.argv[2],Verbose=False)

