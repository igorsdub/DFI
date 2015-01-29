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
if len(sys.argv) < 2: 
    print __doc__
    sys.exit()

def colorbydfi(CSVFIL,PDBFIL,Verbose=False):
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
    pdbid = CSVFIL.split('-')[0]

    data = pd.read_csv(CSVFIL)
    print "Reading in: %s"%(CSVFIL)

    ATOMS = []
    io.pdb_reader(PDBFIL,ATOMS)
    print "Reading in %s"%(PDBFIL)

    print "Adding b-factors"
    for i in range(len(ATOMS)):
        if ATOMS[i].atom_name.strip(' ') == 'CA':
            resind = ATOMS[i].res_index 
            chainind = ATOMS[i].chainID
            val= data[(data.ResI == int(resind)) & (data.ChainID == chainind.strip(' '))].pctdfi.values[0]
            if Verbose:
                print ATOMS[i].res_index, ATOMS[i].temp_factor,val
            ATOMS[i].temp_factor = val

    io.pdb_writer(ATOMS,filename=pdbid+'-dficolor.pdb')






