#!/usr/bin/env python 
"""
Drop Columns Unnamed 0, Unammed 0.1, Unnamed 0.1
==================================================

Usage
-----
dropcols.py CSVFIL

"""

import sys 
if len(sys.argv) < 2:
    print __doc__
    exit(1)

def drop_cols():
    import pandas as pd
    ls_dropcols = ['Unnamed: 0', 'Unnamed: 0.1', 'Unnamed: 0.1']

    fname=sys.argv[1]
    pd.read_csv(fname,index_col='ResI')
    data=pd.read_csv(fname,index_col='ResI')
    for dc in ls_dropcols:
        if dc in data.columns:
            data.drop(dc, axis=1, inplace=True)
    data.to_csv(fname)

if __name__ == "__main__":
    drop_cols()
