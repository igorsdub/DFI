import numpy as np
import dfi
from dfi.datafiles import example_pdb


def test_pctrank_ascending():
    a = np.array([1, 2, 3, 4, 5])
    assert np.all(dfi.dfi_calc.pctrank(a) == np.array(
        [0.2, 0.4, 0.6, 0.8, 1.0], dtype=float))


def test_pctrank_descending():
    a = np.array([1, 2, 3, 4, 5])
    assert np.all(dfi.dfi_calc.pctrank(a, inverse=True) ==
                  np.array([1., 0.8, 0.6, 0.4, 0.2], dtype=float))


def test_comlineargs():
    comline = '--fdfi A19 A10 --covar mwcovar.dat --pdb %s' % example_pdb
    dict_parms = dfi.dfi_calc.CLdict(comline.split())
    assert dict_parms['--pdb'] == example_pdb
    assert dict_parms['--covar'] == 'mwcovar.dat'
    assert np.all(dict_parms['--fdfi'] == np.array(['A19', 'A10']))
