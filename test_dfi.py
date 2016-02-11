

##test pctrank code 
def test_pctrank_ascending():
    import numpy as np 
    import dfi 
    a = np.array([1,2,3,4,5])
    assert np.all( dfi.pctrank(a) == np.array([0.2,0.4,0.6,0.8,1.0],dtype=float) )

def test_pctrank_descending():
    import numpy as np 
    import dfi 
    a = np.array([1,2,3,4,5])
    assert np.all( dfi.pctrank(a,inverse=True) == np.array([1.,0.8,0.6,0.4,0.2],dtype=float) )

def test_comlineargs():
    import dfi 
    import numpy as np 
    comline='./dfi.py --fdfi A19 A10 --hess mwcovar.dat --pdb 1l2y.pdb'
    dict_parms = dfi.parseCommandLine( comline.split() )
    
    assert dict_parms['--pdb'] == '1l2y.pdb'
    assert dict_parms['--hess'] == 'mwcovar.dat'
    assert np.all( dict_parms['--fdfi'] == np.array(['A19', 'A10']) )
