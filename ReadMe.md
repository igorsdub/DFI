# Dynamic Flexibility Index #

The DFI Code is writen in python.  

[ ![Codeship Status for avishekkumar/dfi](https://codeship.com/projects/2216d090-addf-0133-68b3-42dfb775ebd5/status?branch=master)](https://codeship.com/projects/132158)

### What is this repository for? ###

* Repository of for storing all versions of the dfi code. 
* 0.0.1
* [Repository](https://bitbucket.org/avishekkumar/dfi)


### Dependencies ###

* NumPy
* SciPy
* Pandas
* Pytest 

### Description ###

DFI Calculates the dynamics functional index. 
Right now cacluates the hessian and inverts it 
and write out to the file pinv_svd.debug. 

### Usage ###
```
dfi.py --pdb PDBFILE [--hess HESSFILE] [--fdfi RESNUMS] --help   
```
### Input ###

* PDBFILE:     PDBFILE
* RESNUMS:     e.g., "1,5,6,8"
* HESSFILE:    Covariance (Inverse Hessian) Matrix in a [NxN] ascii format 
* RESNUMS:     Chain + Residues number in the pdb, e.g. A15 B21

### Output Files ###

* Structure used for DFI: dfi-dficolor.pdb 
* Invert the Hessian: pinv_svd.debug
* Master DFI: dfianalysis.csv 

### Developers ###
* Avishek Kumar avishek.kumar@asu.edu


### References ###

1. Atilgan AR, Durell SR, Jernigan RL, Demirel MC, Keskin Bahar I, Biophys. J., 80:505-15, 2001 
2. Glembo T.J, M.F Thorpe, D.W. Farrell, Z. N. Gerek, and S.B. Ozkan. Collective Dynamics Differentiates Functional Divergence in Protein Evolution. 
PLos Computational Biology, 2012  
3. Hayward, S. and B.L. de Groot, Normal modes and essential dynamics. Methods in Molecular Biology 443:89-106,2008
4. Bahar I, Atilgan AR, Erman B, Fold. & Des., 2:173-81, 1997 
