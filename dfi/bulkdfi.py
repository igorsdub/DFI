import dfi
import glob
from six.moves import zip

def bulk_dfi(pdbfile, Verbose=False):
    """
    Calculates dfi for all covariance matrices. 
    Globs all covariance matrices with _mwcovarmat.dat
    and replaces with -dfianalysis.csv for the output files. 

    Input
    -----
    pdbfile: fname
       Name of pdbfile to be used in DFI calculations 

    Output
    ------
    -dfianalysis.csv: fname
       Writes out to file 
    """
    covar_files = glob.glob('*_mwcovarmat.dat')
    dfi_files = map(lambda x: x.replace(
        '_mwcovarmat.dat', '-dfianalysis.csv'), covar_files)
    for covar, dfi_fil in zip(covar_files, dfi_files):
        dfi.calc_dfi(pdbfile, covar=covar,
                     writetofile=True, dfianalfile=dfi_fil)
