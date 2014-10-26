#!/usr/bin/env python 
"""
DFI Calculates the dynamics functional index. 
Right now cacluates the hessian and iverts it 
and write out to the file pinv_svd.debug. 

USAGE: dfi.py PDB

PDB:    PDB FILE     

"""

import sys 
import pdbmunge 
import pdbio 
import numpy as np 
from scipy import linalg as LA


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

def calchessian(resnum,x,y,z,gamma=10,Verbose=False):
    """ Calculates the hessian and retuns the result """

    print "Calculating the Hessian..."
    #gamma=10 
    numresthree = 3*resnum 
    hess = np.zeros((numresthree,numresthree))
  
       
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
                print "i:%d j:%f"%(i,j)
                print "x_i:%f y_i:%f z_i:%f"%(x_i,y_i,z_i)
                print "x_j:%f y_j:%f z_j:%f"%(x_j,y_j,z_j)
                print "x_ij:%f y_ij:%f z_ij:%f r:%f sprngcnst:%f"%(x_ij,y_ij,z_ij,r,sprngcnst)
                       

            #creation of Hii 
            hess[3*i,3*i] += sprngcnst*(x_ij*x_ij/r) 
            hess[3*i+1,3*i+1] += sprngcnst*(y_ij*y_ij/r) #(2,2) (5,5)
            hess[3*i+2,3*i+2] += sprngcnst*(z_ij*z_ij/r)     #(3,3) (6,6)

            hess[3*i,3*i+1] += sprngcnst*(x_ij*y_ij/r) #(1,2) (4,5)
            hess[3*i,3*i+2] += sprngcnst*(x_ij*z_ij/r)   #(1,3) (4,6)
            hess[3*i+1,3*i] += sprngcnst*(y_ij*x_ij/r) #(2,1) (5,4)
             
            hess[3*i+1,3*i+2] += sprngcnst*(y_ij*z_ij/r)   #(2,3) (5,6)
            hess[3*i+2,3*i] += sprngcnst*(x_ij*z_ij/r)   #(3,1) (6,4)
            hess[3*i+2,3*i+1] += sprngcnst*(y_ij*z_ij/r)   #(3,2) (6,5)
            
            #creation of Hij 
            hess[3*i,3*j] -= sprngcnst*(x_ij*x_ij/r) #(1,4) (4,1)
            hess[3*i+1,3*j+1] -= sprngcnst*(y_ij*y_ij/r) #(2,5) (5,2)
            hess[3*i+2,3*j+2] -= sprngcnst*(z_ij*z_ij/r)     #(3,6) (6,3)
            
            hess[3*i,3*j+1] -= sprngcnst*(x_ij*y_ij/r) #(1,5) (4,2)
            hess[3*i,3*j+2] -= sprngcnst*(x_ij*z_ij/r)   #(1,6) (4,3)
            hess[3*i+1,3*j] -= sprngcnst*(y_ij*x_ij/r) #(2,4) (5,1)
             
            hess[3*i+1,3*j+2] -= sprngcnst*(y_ij*z_ij/r)   #(2,6) (5,3)
            hess[3*i+2,3*j] -= sprngcnst*(x_ij*z_ij/r)   #(3,4) (6,1)
            hess[3*i+2,3*j+1] -= sprngcnst*(y_ij*z_ij/r)   #(3,5) (6,2)

    print "Finished Calculating the Hessian..."
    return hess  

def flatandwrite(matrix,outfile):
    outfile=open(outfile,'w')
    for f in matrix.flatten():
        outfile.write('%f\n'%f)
    outfile.close()
    


if __name__ == "__main__":
    Verbose = False #Setting for Debugging  
        
    #parse the input 
    #Add a check to make sure it is a file if not then just download from the pdb. 
    pdbid = sys.argv[1] 
    ATOMS = [] 
    pdbio.pdb_reader(pdbid,ATOMS,CAonly=True,chainA=True)
    pdbio.pdb_writer(ATOMS,msg="HEADER dfi target, CAonly and chainA",filename='dfi-out.pdb')
    x,y,z = getcoords(ATOMS) 
    
    #start computing the Hessian 
    numres = len(ATOMS)
    numresthree = 3 * numres
    hess = calchessian(numres,x,y,z,Verbose)
    e_vals, e_vecs = LA.eig(hess)
    if(Verbose):
        flatandwrite(hess,'hesspy.debug')
    
    #print out eigenvalues 
    i=1
    outfile=open('eigenvalues.txt','w')
    for val in np.sort(e_vals):
        outfile.write("%d\t%f\n"%(i,val)) 
        i += 1
    outfile.close()
        

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
    pinv_svd = np.dot(np.dot(U,np.diag(invw)),Vt)
    flatandwrite(pinv_svd,'pinv_svd.debug')
    print "Hessian inverted and written out to pinv_svd.debug"
    

    exit()
    

    
    

    
    
            

