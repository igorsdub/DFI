# Dynamic Flexibility Index #



[ ![Codeship Status for avishekkumar/dfi](https://codeship.com/projects/2216d090-addf-0133-68b3-42dfb775ebd5/status?branch=master)](https://codeship.com/projects/132158)




### Dependencies ###

See the requirements file [dependencies](https://github.com/avishek-r-kumar/dfi/requirements.txt)

### Description ###

DFI Calculates the dynamic functional index. 

### How to install ###
You can clone this repo:
```
git clone https://github.com/avishek-r-kumar/dfi.git
```
OR to get the latest commit 
```
git clone --depth 1 https://github.com/avishek-r-kumar/dfi.git
```
To install the dependencies you can use pip
```
pip install -r requirements.txt 
```

### Adding the DFI directory to the PATH (Bash) ###
You can also import the module from anywhere in your filesystem by appending the 
dfi directory to the PYTHONPATH variable, e.g.,
```
export PYTHONPATH=${PYTHONPATH}:~/dfi #assuming this is where the dfi fold is in your filesystem 
```
then you should be able to `import dfi_calc` or `import ColorDFI` or `import dfi` anywhere
in the filesystem. 
 

### Usage ###
```
dfi_calc.py --pdb PDBFILE [--hess HESSFILE] [--fdfi RESNUMS] --help   
```
OR
```
dfi.py [UNIPROT IDS]
```
### Example ###
#### Run just bare DFI on a protein ####
```
./dfi_calc.py --pdb 1l2y.pdb --fdfi A10 
```
This will run dfi on 1l2y.pdb and write out to 1l2y-dfianalysis.csv,
1l2y-fdficolor.pdb, and 1l2y-dficolor.pdb.
	
#### Run based on UniprotID ####
```
./dfi.py P42771
```
This will take the UniprotID and blast it on the NCBI server, find the
highest PDB hit and then run dfi analysis on that PDB. The output is
P42771-1DC2-dfianalysis.csv. 
*Note: If you query the NCBI server too often it will push your query
down the queue*

### Input ###

* PDBFILE:     PDBFILE
* HESSFILE:    Covariance (Inverse Hessian) Matrix in a [NxN] ascii format 
* RESNUMS:     Chain + Residues number in the pdb, e.g. A15 B21

### Output Files ###

* Structure used for DFI: dfi-dficolor.pdb 
* Master DFI: dfianalysis.csv 

### Developers ###
* Avishek Kumar avishek.kumar@asu.edu


### References ###

1. Atilgan AR, Durell SR, Jernigan RL, Demirel MC, Keskin Bahar I, Biophys. J., 80:505-15, 2001 
2. Glembo T.J, M.F Thorpe, D.W. Farrell, Z. N. Gerek, and S.B. Ozkan. Collective Dynamics Differentiates Functional Divergence in Protein Evolution. 
PLos Computational Biology, 2012  
3. Hayward, S. and B.L. de Groot, Normal modes and essential dynamics. Methods in Molecular Biology 443:89-106,2008
4. Bahar I, Atilgan AR, Erman B, Fold. & Des., 2:173-81, 1997
5. Kumar, A., Glembo, T. J. & Ozkan, S. B. The Role of Conformational Dynamics and Allostery in the Disease Development of Human Ferritin. Biophysical Journal (2015). doi:10.1016/j.bpj.2015.06.060
6. Butler, B. M., Gerek, Z. N., Kumar, S. & Ozkan, S. B. Conformational dynamics of nonsynonymous variants at protein interfaces reveals disease association: The Role of Dynamics in Neutral and Damaging nsSNVs. Proteins: Structure, Function, and Bioinformatics 83, 428–435 (2015).
7. Kumar, A., Butler, B. M., Kumar, S. & Ozkan, S. B. Integration of structural dynamics and molecular evolution via protein interaction networks: a new era in genomic medicine. Current Opinion in Structural Biology 35, 135–142 (2015).

 
