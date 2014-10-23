#!/usr/bin/env python 
"""
USAGE: dfi.py PDBID

PDBID:    Protein DataBank Code 

"""

import sys 
import pdbmunge 
import pdbio 
import numpy as np 


if len(sys.argv) < 2:
    print __doc__ 
    exit()


def getcoords(ATOMS):
    """ Returns x,y and z numpy arrays of coordinates """
    x = []
    y = []
    z = [] 

    for atom in ATOMS:
        if(Verbose):
            print atom.atom_name
            print "%s %f %f %f"%(atom.atom_name,atom.x,atom.y,atom.z)
        if atom.atom_name == 'CA ':
            x.append(atom.x)
            y.append(atom.y)
            z.append(atom.z)
    
    x = np.array(x,dtype=float)
    y = np.array(y,dtype=float)
    z = np.array(z,dtype=float)
    
    return x,y,z 

def calchessian(hess,x,y,z):
    """ Calculates the hessian and retuns the result """
    return 0 


if __name__ == "__main__":
    Verbose = True 
    gamma = 100
    

    #parse the input 
    pdbid = sys.argv[1] 
    
    #get the pdb, then extract the chain and then the alpha carbons and return the name of the file 
    fname = pdbmunge.extractA_CA(pdbid)
    ATOMS = [] 
    pdbio.pdb_reader(fname,ATOMS)

    x,y,z = getcoords(ATOMS) 

    #start computing the Hessian 
    numres = len(ATOMS)
    numresthree = 3 * numres
    

    hess = np.zeros((numresthree,numresthree))
    

    print hess 
    
    #compute the Hessian 
    #compute the Hii terms 
    for i in len(numres):
        for j in len(numres):
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
            sprngcnst = (gamma*gamma*gamm)/(r*r*r)
                       
            #creation of Hii 
            hess[3*i-2,3*i-2] += sprngcnst*(x_ij*x_ij/r)
            hess[3*i-1,3*i-1] += sprngcnst*(y_ij*y_ij/r)
            hess[3*i,3*i] += sprngcnst*(z_ij*z_ij/r)

            hess[3*i-2,3*i-1] += sprngcnst*(x_ij*y_ij/r)
            hess[3*i-2,3*i] += sprngcnst*(x_ij*z_ij/r)
            hess[3*i-1,3*i-2] += sprngcnst*(y_ij*x_ij/r)
            
            hess[3*i-1,3*i] += sprngcnst*(y_ij*z_ij/r)
            hess[3*i,3*i-2] += sprngcnst*(x_ij*z_ij/r)
            hess[3*i,3*i-1] += sprngcnst*(y_ij*z_ij/r)
            
            #creation of Hij 
            hess[3*i-2,3*j-2] -= sprngcnst*(x_ij*x_ij/r)
            hess[3*i-1,3*j-1] -= sprngcnst*(y_ij*y_ij/r)
            hess[3*i,3*j] -= sprngcnst*(z_ij*z_ij/r)
            
            hess[3*i-2,3*j-1] -= sprngcnst*(x_ij*y_ij/r)
            hess[3*i-2,3*j] -= sprngcnst*(x_ij*z_ij/r)
            hess[3*i-1,3*j-2] -= sprngcnst*(y_ij*x_ij/r)
            
            hess[3*i-1,3*j] -= sprngcnst*(y_ij*z_ij/r)
            hess[3*i,3*j-2] -= sprngcnst*(x_ij*z_ij/r)
            hess[3*i,3*j-1] -= sprngcnst*(y_ij*z_ij/r)

            #need to check is the Hessian is going to work 
    
            

