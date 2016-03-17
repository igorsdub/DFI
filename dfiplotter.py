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


def plotdfi(df, quant, pdbid):
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd

    sns.set_style("whitegrid")
    plt.figure(figsize=(22, 12))
    sns.set_context("poster", font_scale=1.25, rc={
                    "lines.linewidth": 1.25, "lines.markersize": 8})

    df[quant].plot(marker='o', label='pdbid')

    plt.ylabel('%DFI')
    plt.xlabel('Residue Index')
    plt.savefig(pdbid + '-' + quant + '.png')
    return plt
    print "Printed: %s" % (pdbid + '-' + quant + '.png')
