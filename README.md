# Dynamic Flexibility Index #


### Short Description ###
DFI is a python package that calculates protein dynamics from the
[Ozkan Group][OzkanLab] at Arizona State University.

[OzkanLab]: <http://ozkanlab.physics.asu.edu> "Ozkan Lab Website"

### Longer Description ###
The Dynamic Flexiblity index (DFI) is a measure of each residue's contribution to
a protein's dynamics. A low %DFI score indicates a rigid portion of the protein
while a high %DFI score indicates a flexible portion of the protein. The %DFI
score has been used as a predictive feature in genetic disease prediction. 
Mutations in genomes can lead to proteins which can misfunction, manifesting in 
genetic disease. Our proteome-wide analysis using DFI indicates that certain 
sites play a critical role in the functinally reated dynamics (i.e, those with
low dfi values); therefore, mutations at those sites are more likely to be 
associated with disease. 

### Dependencies ###

See the requirements file [dependencies][Requirements]

[Requirements]: <https://raw.githubusercontent.com/avishek-r-kumar/DFI/master/requirements.txt>

## How to install ##
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

## Usage ##
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

1. Atilgan AR, Durell SR, Jernigan RL, Demirel MC, Keskin Bahar I, Biophys. J., 
80:505-15, 2001 
2. Glembo T.J, M.F Thorpe, D.W. Farrell, Z. N. Gerek, and S.B. Ozkan. Collective
 Dynamics Differentiates Functional Divergence in Protein Evolution. 
PLos Computational Biology, 2012  
3. Hayward, S. and B.L. de Groot, Normal modes and essential dynamics. Methods i
n Molecular Biology 443:89-106,2008
4. Bahar I, Atilgan AR, Erman B, Fold. & Des., 2:173-81, 1997
5. Kumar, A., Glembo, T. J. & Ozkan, S. B. The Role of Conformational Dynamics a
nd Allostery in the Disease Development of Human Ferritin. Biophysical Journal (
2015). doi:10.1016/j.bpj.2015.06.060
6. Butler, B. M., Gerek, Z. N., Kumar, S. & Ozkan, S. B. Conformational dynamics
 of nonsynonymous variants at protein interfaces reveals disease association: Th
e Role of Dynamics in Neutral and Damaging nsSNVs. Proteins: Structure, Function
, and Bioinformatics 83, 428–435 (2015).
7. Kumar, A., Butler, B. M., Kumar, S. & Ozkan, S. B. Integration of structural 
dynamics and molecular evolution via protein interaction networks: a new era in 
genomic medicine. Current Opinion in Structural Biology 35, 135–142 (2015).

 
