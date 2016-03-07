#!/usr/bin/env python 
"""
DFI (Dynamic Flexibility Index)
===============================

Description
------------
DFI Calculates the dynamics flexibility index. 
Program calculates the hessian and inverts it 
and write out to the file pinv_svd.debug. Does
a-dfi as well now.

Requirements
------------
Python 2.7.5
NumPy 1.4
SciPy 1.4 
Pandas 

Usage
-----
dfi.py --pdb PDBFILE [--hess HESSFILE] [--chain CHAINID] [--fdfi RESNUMS] --help   

Input
-----
PDBFILE:     PDBFILE
RESNUMS:     e.g., "1,5,6,8"
HESSFILE:    Covariance (Inverse Hessian) Matrix in a [NxN] ascii format  
RESNUMS:     Chain + Residues number in the pdb, e.g. A15 B21

Output 
------
* Structure used for DFI: dfi-out.pdb 
* Inverted Hessian: pinv_svd.debug 
* Master DFI: dfianalysis.csv      

Example
-------
```
./dfi.py --pdb 1a2x_BA_1.pdb --hess covariance.dat --chain A --fdfi A15 A95 A98 A101 A102 A118 A119 A126 B17 B20 B21 B22 B24 B29
```
"""

import sys 
import pdbio 
import os 
import numpy as np 
import pandas as pd 
import ColorDFI
import dfiplotter

from scipy import linalg as LA
from scipy import stats 


def getcoords(ATOMS,Verbose=False):
    """ Returns x,y and z numpy arrays of coordinates """
    x = []
    y = []
    z = [] 
    bfac = []

    for atom in ATOMS:
        if(Verbose):
            print atom.atom_name
            print "%s %f %f %f %f"%(atom.atom_name,atom.x,atom.y,atom.z,atom.temp_factor)
        if atom.atom_name == 'CA ':
            x.append(atom.x)
            y.append(atom.y)
            z.append(atom.z)
            bfac.append(atom.temp_factor)

    x = np.array(x,dtype=float)
    y = np.array(y,dtype=float)
    z = np.array(z,dtype=float)
    bfac = np.array(bfac,dtype=float)

    return x,y,z,bfac  

def calchessian(resnum,x,y,z,gamma=float(100),cutoff=None,Verbose=False):
    """ 
    Calculates the hessian and retuns the result 
    ============================================
    
    calchessian(resnum,x,y,z,gamma=100,cutoff=None,Verbose=False)

    Inputs 
    ------
    resnum: number of residues 
    x: array of x coordinates
    y: array of y coordinates
    z: array of z coordinates
    gamma: value of spring constant(default set to 100)
    cutoff: value of cutoff when using a distance based Hessian (default None)
    Verbose: Verbose Output for debug mode (default False).
    
    Output 
    ------
    hess: numpy array of the Hessian 
    """

    print "Calculating the Hessian..."
    
    numresthree = 3*resnum 
    hess = np.zeros((numresthree,numresthree))
    gamma=100 
       
    #compute the Hessian 
    #compute the Hii terms 
    for i in range(resnum):
        for j in range(resnum):
            if i==j:
                continue
            x_i = x[i]
            y_i = y[i]
            z_i = z[i]
            x_j = x[j]
            y_j = y[j]
            z_j = z[j]
            x_ij = x_i - x_j 
            y_ij = y_i - y_j 
            z_ij = z_i - z_j 
            r = x_ij*x_ij + y_ij*y_ij + z_ij*z_ij 
            sprngcnst = (gamma*gamma*gamma)/(r*r*r)
            if(Verbose):
                print "sprngcnst:",sprngcnst,"gamma",gamma
            if(cutoff):
                if sqrt(r) > cutoff:
                    if(Verbose):
                        print cutoff
                    #sprngcnst = 0. 
            if(Verbose):
                print "i:%d j:%f"%(i,j)
                print "x_i:%f y_i:%f z_i:%f"%(x_i,y_i,z_i)
                print "x_j:%f y_j:%f z_j:%f"%(x_j,y_j,z_j)
                print "x_ij:%f y_ij:%f z_ij:%f r:%f sprngcnst:%f"%(x_ij,y_ij,z_ij,r,sprngcnst)
                       

            #creation of Hii 
            hess[3*i,3*i] += sprngcnst*(x_ij*x_ij/r) 
            hess[3*i+1,3*i+1] += sprngcnst*(y_ij*y_ij/r) 
            hess[3*i+2,3*i+2] += sprngcnst*(z_ij*z_ij/r) 

            hess[3*i,3*i+1] += sprngcnst*(x_ij*y_ij/r) 
            hess[3*i,3*i+2] += sprngcnst*(x_ij*z_ij/r) 
            hess[3*i+1,3*i] += sprngcnst*(y_ij*x_ij/r) 
             
            hess[3*i+1,3*i+2] += sprngcnst*(y_ij*z_ij/r) 
            hess[3*i+2,3*i] += sprngcnst*(x_ij*z_ij/r)   
            hess[3*i+2,3*i+1] += sprngcnst*(y_ij*z_ij/r) 
            
            #creation of Hij 
            hess[3*i,3*j] -= sprngcnst*(x_ij*x_ij/r) 
            hess[3*i+1,3*j+1] -= sprngcnst*(y_ij*y_ij/r) 
            hess[3*i+2,3*j+2] -= sprngcnst*(z_ij*z_ij/r) 
            
            hess[3*i,3*j+1] -= sprngcnst*(x_ij*y_ij/r) 
            hess[3*i,3*j+2] -= sprngcnst*(x_ij*z_ij/r) 
            hess[3*i+1,3*j] -= sprngcnst*(y_ij*x_ij/r) 
             
            hess[3*i+1,3*j+2] -= sprngcnst*(y_ij*z_ij/r) 
            hess[3*i+2,3*j] -= sprngcnst*(x_ij*z_ij/r)   
            hess[3*i+2,3*j+1] -= sprngcnst*(y_ij*z_ij/r) 

    print "Finished Calculating the Hessian..."
    return hess  

def flatandwrite(matrix,outfile):
    """Flattens out a matrix to a Nx1 column and write out to a file. """
    outfile=open(outfile,'w')
    for f in matrix.flatten():
        outfile.write('%f\n'%f)
    outfile.close()
    
def dfianal(fname,Array=False):
    """Calculate various dfi quantities and then output"""
    if not(Array):
        with open(fname,'r') as infile:
            dfi = np.array([x.strip('\n') for x in infile], dtype=float )
    else:
        dfi = fname 
    dfirel = dfi/np.mean(dfi)
    dfizscore = stats.zscore(dfi)
    dfiperc = [] 
    lendfi = float(len(dfi))
    for m in dfi:
        amt = np.sum(dfi <=m)
        dfiperc.append(amt/lendfi)
    dfiperc = np.array(dfiperc,dtype=float)
    
    return dfi, dfirel, dfiperc, dfizscore 

def pctrank(dfi,inverse=False):

    if type(dfi).__module__ != 'numpy':
        raise ValueError('Input needs to be a numpy array')

    dfiperc = [] 
    lendfi = float(len(dfi))
    
    for m in dfi:
        
        if inverse:
            amt = np.sum(dfi >=m)
        else:
            amt = np.sum(dfi <=m)
        dfiperc.append(amt/lendfi)

    return np.array(dfiperc,dtype=float)

def calcperturbMat(invHrs,direct,resnum,Normalize=True):
    perturbMat = np.zeros((resnum,resnum))
    for k in range(len(direct)):
        peturbDir = direct[k,:]
        for j in range(int(resnum)):
            delforce = np.zeros(3*resnum)
            delforce[3*j:3*j+3] = peturbDir 
            delXperbVex = np.dot(invHrs,delforce)
            delXperbMat = delXperbVex.reshape((resnum,3))
            delRperbVec = np.sqrt(np.sum(delXperbMat*delXperbMat,axis=1))
            perturbMat[:,j] += delRperbVec[:]
    perturbMat /= 7

    if(Normalize):
        nrmlperturbMat = perturbMat/np.sum(perturbMat)
    else:
        print "WARNING: The perturbation matrix is not NORMALIZED"
        nrmlperturbMat = perturbMat
    
    return nrmlperturbMat 


def CLdict(argv):
    """Returns a dictionary of the command lines options from sys.argv"""
    comline_arg={}
    for s in argv:
        if s ==  "--pdb":
            ind = argv.index(s)
            comline_arg[s] = argv[ind+1]
            if (os.path.isfile(argv[ind+1]) != True):
                print "File "+ argv[ind+1] +" not found."
                sys.exit(1)

        if s ==  "--hess":
            ind = argv.index(s)
            comline_arg[s] = argv[ind+1]
                    
        if s == "--chain":
            ind = argv.index(s)
            comline_arg[s] = argv[ind+1]
            
               
        if s ==  "--fdfi":
            ind = argv.index(s)
            resvals = []
            for res in argv[ind+1:]:
                if res.startswith("--"):
                    break 
                else:
                    resvals.append(res)
            comline_arg[s] = np.array(resvals)
            print resvals 

        if s == "--help":
            print __doc__
            sys.exit(1)
            
    if ("--pdb" not in argv):
        print argv
        print "No --pdb"
        print __doc__
        sys.exit(1)
            
    return comline_arg

def chainresmap(ATOMS,Verbose=False):
    """
    Returns a dict object with the chainResNum as the key and the index
    of the atom 
    """
    table = {}
    for i in range(len(ATOMS)):
        if ATOMS[i].chainID==' ':
            entry = ATOMS[i].res_index
        else:
            entry = ATOMS[i].chainID + ATOMS[i].res_index #str(ATOMS[i].res_index)
        table[entry] = i
    if(Verbose):
        print table 
    return table 

def fdfiresf(ls_chain,table):
    """Returns numpy array of f-dfi res"""
    ls_ind = []
    
    for res in ls_chain:
        if res in table:
            ls_ind.append( table[res] )
        else:
            print "WARNING: Can't find %s"%res 
            continue 
    return np.array(ls_ind,dtype=int)


def fdfires_cords(fdfires,x,y,z):
    """
    fdfires_cords(fdfires,x,y,z)
    fdfires: indices of fdfires 
    x: array of x-coordinates 
    y: array of y-coordinates 
    z: array of z-coordinates
    return: nx3 matrix of f-DFI coordinates 
    """
    print x[:10],y[:10],z[:10]
    return np.column_stack((x[fdfires],y[fdfires],z[fdfires]))

def rdist(r,fr):
    """
    rdist
    r = array of coordinates 
    fr = nx3 matrix of f-DFI coordinates 
    return rdist: array of distances from f-DFI sites 
    """
    r_ij = fr - r
    rr = r_ij*r_ij
    #print "r_ij",r_ij,"rr",rr,"r",np.sqrt(rr.sum(axis=1))
    return np.sqrt(rr.sum(axis=1))

def outputToDF(ATOMS,dfi,pctdfi,fdfi=None,pctfdfi=None,ls_ravg=None,ls_rmin=None,outfile=None):
   """
   Outputs the results of the DFI calculation to a DataFrame and csv file 
   """
   mapres={'ALA':'A',
            'CYS':'C',
            'ASP':'D',
            'GLU':'E',
            'PHE':'F',
            'GLY':'G',
            'HIS':'H',
            'ILE':'I',
            'LYS':'K',
            'LEU':'L',
            'MET':'M',
            'PRO':'P',
            'ARG':'R',
            'GLN':'Q',
            'ASN':'N',
            'SER':'S',
            'THR':'T',
            'TRP':'W',
            'TYR':'Y',
            'VAL':'V'}
   dfx = pd.DataFrame()
   dfx['ResI'] = [ ATOMS[i].res_index.strip(' ') for i in xrange(len(ATOMS))]    
   dfx['dfi'] = dfi 
   dfx['pctdfi'] = pctdfi 
   dfx['ChainID'] = [ATOMS[i].chainID for i in xrange(len(ATOMS))]
   dfx['Res'] = [ATOMS[i].res_name for i in xrange(len(ATOMS))]
   dfx['R'] = dfx['Res'].map(mapres)


   if type(fdfi).__module__ == 'numpy':
       dfx['fdfi'] = fdfi 
       dfx['pctfdfi'] = pctfdfi 
       dfx['ravg'] = ls_ravg
       dfx['rmin'] = ls_rmin 
       mask = (dfx['rmin'] > 8.0) & (dfx['pctfdfi'] > 0.75)
       dfx['A'] = mask.map(lambda x: 'A' if x else 'NotA')
        
   if(outfile):
       dfx.to_csv(outfile,index=False)
   
   return dfx 

def top_quartile_pos(pctfdfi,rlist):
    """
    returns a list of indices of positions in the top quartile of pctdfi 
    """
    return [i for i,val in enumerate(pctfdfi) if val > 0.75]        
        

if __name__ == "__main__" and len(sys.argv) < 2:
    print __doc__
    sys.exit(1)

def parseCommandLine(argv):
    """
    Parse command lines input
    
    Input 
    -----
    ls of command line input 

    Output
    ------
    pdbfile: name of pdb file to run dfi calculation 
    pdbid: 4 Letter PDB code 
    mdhess: name of file that contains hessian matrix from MD 
    ls_reschain: list of f-dfi Residues (e.g., [A17,A19])
    chain_name: list of chain to select (Depracated)
    
    """
    comlinargs=CLdict(argv)
    pdbfile = comlinargs['--pdb']
    pdbid = pdbfile.split('.')[0]
    mdhess=bool( comlinargs.get('--hess',"") )
    ls_reschain=comlinargs.get('--fdfi',[])
    chain_name = comlinargs.get('--chain','A')
    
    return pdbfile, pdbid, mdhess, ls_reschain, chain_name


def dfi(argv):
    Verbose = False #Setting for Debugging  
    #Parse the input 
   
    comlinargs=CLdict(argv)
      
    pdbfile = comlinargs['--pdb']
    pdbid = pdbfile.split('.')[0]
    mdhess=bool( comlinargs.get('--hess',"") )
    ls_reschain=comlinargs.get('--fdfi',[])
    chain_name = comlinargs.get('--chain','A')

    #output file name 
    eigenfile = pdbid+'-eigenvalues.txt'
    invhessfile = pdbid+'-pinv_svd.debug'
    dfianalfile = pdbid+'-dfianalysis.csv'
     

    #read in the pdb file 
    ATOMS = [] 
    pdbio.pdb_reader(pdbfile,ATOMS,CAonly=True,noalc=True,chainA=False,chain_name=chain_name,Verbose=False)
    x,y,z,bfac = getcoords(ATOMS) 

    #parse the f-dfi inputs 
    
    if len(ls_reschain) > 0:
        print "f-dfires"
        fdfiset = set(ls_reschain)
        ls_reschain = list(fdfiset)
        ls_reschain.sort()   
        print ls_reschain
        print "Number of f-dfi: %d" %len(ls_reschain)
        fdfires = np.sort( fdfiresf(ls_reschain,chainresmap(ATOMS)) )
    else:
        fdfires = np.array([],dtype=int) 

    print "fdfires numpy"
    print fdfires

    
    #start computing the Hessian 
    numres = len(ATOMS)
    numresthree = 3 * numres 

    if not(mdhess):
        #numres = len(ATOMS)
        #numresthree = 3 * numres
        hess = calchessian(numres,x,y,z,Verbose)
        e_vals, e_vecs = LA.eig(hess)
        if(Verbose):
            print "Hessian"
            print hess
            flatandwrite(hess,'hesspy.debug')
    
        i=1
        with open(eigenfile,'w') as outfile:
            for val in np.sort(e_vals):
                outfile.write("%d\t%f\n"%(i,np.real(val) ) )
                i += 1
        
        U, w, Vt = LA.svd(hess,full_matrices=False)
        if(Verbose):
            print U.shape
            print w.shape 
            print Vt.shape 

        S = LA.diagsvd(w,len(w),len(w))
        print "Checking If the SVD went well..."
        print np.allclose(hess,np.dot(U,np.dot(S,Vt)))
     
        if(Verbose):
            flatandwrite(U,'Upy-test.debug')
            flatandwrite(w,'wpy-test.debug')
            flatandwrite(Vt,'Vtpy-test.debug')

        #the near zero eigenvalues blowup the inversion so 
        #we will truncate them and add a small amount of bias 
        print "Inverting the Hessian..."
        tol = 1e-6 
        singular = w < tol 
        invw = 1/w
        invw[singular] = 0.
        invHrs = np.dot(np.dot(U,np.diag(invw)),Vt)
        flatandwrite(invHrs,invhessfile)
        print "Number of near-singular eigenvalues: %f"%np.sum(singular)
        print "Hessian inverted and written out to pinv_svd.debug"


    #import the inverse Hessian and turn into a function  
    if not(mdhess):
        print "Reading theinverse Hessian from %s"%(invhessfile)
        with open(invhessfile,'r') as infile:
            invH = np.array([ x.strip('\n') for x in infile],dtype=float) 
       
        resnumsq = len(invH)
        resnum = np.sqrt(resnumsq)/3
        invHrs = invH.reshape((3*resnum,3*resnum),order='F')
    else:
        invHrs=np.loadtxt( comlinargs['--hess'] )
        print "From MD invhess"
        print invHrs
        print invHrs.shape 
    if(Verbose):
        print "invHrs"
        print invHrs 
    
    #RUN DFI CODE HERE 
    print "Creating the perturbation directions"
    directions = np.vstack(([1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1],[1,1,1]))
    normL = np.linalg.norm(directions,axis=1)
    direct=directions/normL[:,None]
    print "Perturbation Directions"
    print direct 


    print "Calculating the peturbation matrix"
    nrmlperturbMat = calcperturbMat(invHrs,direct,numres)
    dfi = np.sum(nrmlperturbMat,axis=1)
    mdfi = np.sum(nrmlperturbMat,axis=0)
        
    
    dfi, reldfi, pctdfi, zscoredfi = dfianal(dfi,Array=True)
    mdfi, relmdfi, pctmdfi, zscoremdfi = dfianal(mdfi,Array=True)

    
    #f-dfi 
    print "Amount of f-dfi res:"+str(len(fdfires))
    print fdfires 
    fdfifile='fdfi-Avg.dat'
    if len(fdfires) > 0:
        #Very clunky impplementtation 
        fdfitop=np.sum(nrmlperturbMat[:,fdfires],axis=1)/len(fdfires)
        fdfibot=np.sum(nrmlperturbMat,axis=1)/len(nrmlperturbMat)
        flatandwrite(fdfitop/fdfibot,fdfifile)
        fdfi,relfdfi,pctfdfi,zscorefdfi = dfianal(fdfifile)

        x,y,z,bfac = getcoords(ATOMS)
        rlist = np.column_stack((x,y,z)) #dump into a list 
        fr = fdfires_cords(fdfires,x,y,z) #get coordinates of the f-dfi residues 
        ls_ravg = np.array([ rdist(r,fr).mean() for r in rlist]) 
	ls_rmin = np.array([ rdist(r,fr).min() for r in rlist])

    if len(fdfires) > 0:
        df_dfi = outputToDF(ATOMS,dfi,pctdfi,fdfi=fdfi,pctfdfi=pctfdfi,ls_ravg=ls_ravg,ls_rmin=ls_rmin,outfile=dfianalfile)
    else:
        df_dfi = outputToDF(ATOMS,dfi,pctdfi,outfile=dfianalfile)

    ColorDFI.colorbydfi(dfianalfile,pdbfile,colorbyparam='pctdfi',outfile=pdbid+'-dficolor.pdb')
    if len(fdfires) > 0:
        ColorDFI.colorbydfi(dfianalfile,pdbfile,colorbyparam='pctfdfi',outfile=pdbid+'-fdficolor.pdb')

    return pdbid,df_dfi 
    

if __name__ == "__main__":
    print sys.argv
    pdbid , df_dfi = dfi(sys.argv)
    dfiplotter.plotdfi(df_dfi,'pctdfi',pdbid)
