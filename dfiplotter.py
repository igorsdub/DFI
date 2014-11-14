#!/usr/bin/env python  
"""
======================
DFI Plots of pctdfi(s)
======================

Usage
-----

./dfiplotter.py pdbid-dfianalysis.csv 

Description
------------

Plots the pctdfis and fulldfi and prints them out.
Nothing fancy with the plots maybe add seaborn. 

Output
------

pdbid-pctdfi.png
pdbid-mpctdfi.png
pdbid-pcthmdfi.png 
pdbid-DFIfig.png 
"""

def plotdfi(data,quant,pdbid):
    plt.plot(data[quant],'bo-')
    plt.title(quant.replace('pct','%'))
    plt.xlabel('Residue Index')
    plt.savefig(pdbid+'-'+quant+'.png')
    plt.close()
    print "Printed: %s"%(pdbid+'-'+quant+'.png')

import sys 
if __name__=="__main__" and len(sys.argv) < 2:
    print __doc__
else:
    
    import numpy as np 
    import matplotlib.pyplot as plt 
    import pandas as pd 
    
    dfifile = sys.argv[1]
    pdbid = dfifile.split('-')[0]

    data = pd.read_csv(dfifile,index_col=0)
    d=data.loc[:,['pctdfi','pctmdfi','pcthmdfi']]
    
    fig, axes = plt.subplots(nrows=3, ncols=1,sharex=True)
    for i, c in enumerate(d.columns):
        d[c].plot( ax=axes[i], figsize=(14, 14), title=c, grid=False,fontsize='xx-small')
    plt.setp([a.get_xticklabels() for a in axes], visible=False)
    plt.savefig(pdbid+'-DFIfig.png', bbox_inches='tight')
    print "Printed: %s"%(pdbid+'DFIfig.png')
    plt.close()

    for f in d.columns:
        plotdfi(d,f,pdbid)
    



