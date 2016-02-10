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

def plotdfi(df,quant,pdbid):
    import numpy as np 
    import matplotlib.pyplot as plt
    import seaborn as sns 
    import pandas as pd 

    sns.set_style("whitegrid")
    plt.figure(figsize=(22, 12))
    sns.set_context("poster", font_scale=1.25, rc={"lines.linewidth": 1.25,"lines.markersize":8})


    df[quant].plot(marker='o',label='pdbid')

    plt.ylabel('%DFI')
    plt.xlabel('Residue Index')
    plt.savefig(pdbid+'-'+quant+'.png')
    return plt 
    print "Printed: %s"%(pdbid+'-'+quant+'.png')

#import sys 
#if __name__=="__main__" and len(sys.argv) < 2:
#    print __doc__
#else:
    
#    import numpy as np 
#    import matplotlib.pyplot as plt 
#    import pandas as pd 
    
#    dfifile = sys.argv[1]
#    pdbid = dfifile.split('-')[0]

#    data = pd.read_csv(dfifile,index_col=0)
#    if 'pctfdfi' in data.columns:
#        d=data.loc[:,['pctdfi','pctmdfi','pcthmdfi','pctfdfi']]
#        colors = ['red','cyan','green','blue']
#        nrows=4
#    else:
#        d=data.loc[:,['pctdfi','pctmdfi','pcthmdfi']]
#        colors = ['red','cyan','green']
#        nrows=3

#    plt.style.use('bmh')
#    fig, axes = plt.subplots(nrows=nrows, ncols=1,sharex=True)
    
#    for i, c in enumerate(d.columns):
#        d[c].plot( ax=axes[i], figsize=(14, 14), marker='o', color=colors[i], 
#                   title=c.replace('pct','%'), grid=False,fontsize='xx-small')
#    plt.setp([a.get_xticklabels() for a in axes], visible=False)
#    plt.savefig(pdbid+'-DFIfig.png', bbox_inches='tight')
#    print "Printed: %s"%(pdbid+'-DFIfig.png')
#    plt.close()

#    for f in d.columns:
#        plotdfi(d,f,pdbid)
    
