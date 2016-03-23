import numpy as np


import os.path, sys
sys.path.append(os.path.join(os.path.dirname(
os.path.realpath(__file__)), os.pardir))


import dfi_calc



# test pctrank code


def test_pctrank_ascending():
    a = np.array([1, 2, 3, 4, 5])
    assert np.all(dfi_calc.pctrank(a) == np.array(
        [0.2, 0.4, 0.6, 0.8, 1.0], dtype=float))


def test_pctrank_descending():
    a = np.array([1, 2, 3, 4, 5])
    assert np.all(dfi_calc.pctrank(a, inverse=True) ==
                  np.array([1., 0.8, 0.6, 0.4, 0.2], dtype=float))


def test_comlineargs():
    comline = '--fdfi A19 A10 --hess mwcovar.dat --pdb ./data/1l2y.pdb'
    dict_parms = dfi_calc.CLdict(comline.split())
    assert dict_parms['--pdb'] == './data/1l2y.pdb'
    assert dict_parms['--hess'] == 'mwcovar.dat'
    assert np.all(dict_parms['--fdfi'] == np.array(['A19', 'A10']))
