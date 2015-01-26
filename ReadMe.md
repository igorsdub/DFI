# Dynamic Flexibility Index #

The DFI Code is writen in python.  

### What is this repository for? ###

* Repository of for storing all versions of the dfi code. 
* 0.0.1
* [Repository](https://bitbucket.org/avishekkumar/dfi)

### Dependencies ###

* NumPy
* SciPy
* Test Cases Coming Soon

### Description ###

DFI Calculates the dynamics functional index. 
Right now cacluates the hessian and inverts it 
and write out to the file pinv_svd.debug. 

### Usage ###

dfi.py --pdb PDBFILE [--hess HESSFILE] [--fdfi RESNUMS] --help   

### Input ###

* PDBFILE:     PDBFILE
* RESNUMS:     e.g., "1,5,6,8"
* HESSFILE:    Covariance (Inverse Hessian) Matrix in a [NxN] ascii format 
* RESNUMS:     Chain + Residues number in the pdb, e.g. A15 B21

### Output Files ###

* Structure used for DFI: dfi-out.pdb 
* Eigenvalues: eigenvalues.txt 
* Invert the Hessian: pinv_svd.debug
* DFI: S1-Avg.dat 
* MDFI: S2-Avg.dat 
* Master DFI: dfianalysis.csv 


### TODO ###

* ~~Create and xml file for all inputs to make parameters more flexible~~ 
* ~~Fix bugs with the alc~~ 
* ~~Change name of output files to reflect pdbname and quantities and log files.~~ 
* ~~Fixed Residus Numbers~~ 
* ~~Create a feature to make dfiplots~~
    ~~* Try making publication quality plots instead of bare matplotlib~~  
* ~~Create a feature to just take in the Hessian from MD~~
    * Improve on the implementation of the MD Hessian 
* ~~Check to see how many near singular eigenvalues there are.~~ 
* ~~Do a check to see if the residue numbers are unique, otherwise renumber them.~~
    * ~~Going to trust the PDB file numbers for now.~~  
* Better documentation of the functions
* Add test cases 
* ~~Allow for ChainID~~
* ~~Allowed for identifying chains in the f-dfi~~
* ~~Fix problems with change ids~~
     * ~~Make sure chain id does not have to be a letter~~ 
     * ~~Allow for chain ids not to be found, eliminate from dictionary~~ 

### Who do I talk to? ###

* Avishek Kumar avishek.kumar@asu.edu


### References ###

1. Atilgan AR, Durell SR, Jernigan RL, Demirel MC, Keskin Bahar I, Biophys. J., 80:505-15, 2001 
2. Glembo T.J, M.F Thorpe, D.W. Farrell, Z. N. Gerek, and S.B. Ozkan. Collective Dynamics Differentiates Functional Divergence in Protein Evolution. 
PLos Computational Biology, 2012  
3. Hayward, S. and B.L. de Groot, Normal modes and essential dynamics. Methods in Molecular Biology 443:89-106,2008
4. Bahar I, Atilgan AR, Erman B, Fold. & Des., 2:173-81, 1997 
