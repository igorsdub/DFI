# IPython log file


def test_hessian():
    
    import os.path, sys
    import dfi
    import numpy as np

    a = np.array([64, 0, 0, -64, 0, 0], dtype=float)
    b = np.array([0, 0, 0, 0, 0, 0], dtype=float)
    c = np.array([-64, 0, 0, 64, 0, 0], dtype=float)

    test_hess = np.vstack((a, b, b, c, b, b))

    x = np.array([0, 5], dtype=float)
    y = np.array([0, 0], dtype=float)
    z = np.array([0, 0], dtype=float)
    resnum = 2
    hess = dfi.calchessian(resnum, x, y, z)
    assert np.all(test_hess == hess)
