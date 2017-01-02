# Dynamic Flexibility Index


[![Build Status](https://travis-ci.com/avishekrk/DFI.svg?token=qr1WKDpoEiNDipEKFzrb&branch=master)](https://travis-ci.com/avishekrk/DFI)
[![codecov](https://codecov.io/gh/avishekrk/DFI/branch/master/graph/badge.svg)](https://codecov.io/gh/avishekrk/DFI)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)



## Short Description
DFI is a python package that calculates protein dynamics from the
[Ozkan Group][OzkanLab] at Arizona State University.

[OzkanLab]: <http://ozkanlab.physics.asu.edu> "Ozkan Lab Website"

## Longer Description
The Dynamic Flexiblity index (DFI) is a measure of each residue's contribution to
a protein's dynamics. A low %DFI score indicates a rigid portion of the protein
while a high %DFI score indicates a flexible portion of the protein. The %DFI
score has been used as a predictive feature in genetic disease prediction.
Mutations in genomes can lead to proteins which can misfunction, manifesting in
genetic disease. Our proteome-wide analysis using DFI indicates that certain
sites play a critical role in the functinally reated dynamics (i.e, those with
low dfi values); therefore, mutations at those sites are more likely to be
associated with disease. DFI has been used to study docking, protein evolution,
and disease evolution (see References).

---

# How to install
DFI can be pip installed with the following command
```
pip install git+http://github.com/avishekrk/dfi.git
```
OR to get the latest commit
```
git clone --depth 1 https://github.com/avishek-r-kumar/dfi.git
```

## Dependencies


DFI was written and tested using Python 2.7 and Pyton 3.5.
See the requirements file for dependencies
[dependencies][Requirements].

[Requirements]: <https://raw.githubusercontent.com/avishek-r-kumar/DFI/master/requirements.txt>

To install the dependencies you can use pip
```
pip install -r requirements.txt
```
---


# Usage
```
dfi_calc.py --pdb PDBFILE [--covar COVARFILE] [--fdfi RESNUMS] [--help]
```
OR
```
dfi.py [UNIPROT IDS]
```
## Example
### Run just bare DFI on a protein
```
./dfi_calc.py --pdb 1l2y.pdb --fdfi A10
```
This will run dfi on 1l2y.pdb and write out to 1l2y-dfianalysis.csv,
1l2y-fdficolor.pdb, and 1l2y-dficolor.pdb.

### Run based on UniprotID
```
./dfi.py P42771
```
This will take the UniprotID and blast it on the NCBI server, find the
highest PDB hit and then run dfi analysis on that PDB. The output is
P42771-1DC2-dfianalysis.csv.
*Note: If you query the NCBI server too often it will push your query
down the queue*

### Tutorial
If you would like to call DFI as a package within python. This
[notebook][TutorialLink] provides a short tutoral.

[TutorialLink]: <https://github.com/avishekrk/DFI/blob/master/tutorial/DFI_Tutorial.ipynb> "Tutorial Link"

---

# Input and Output files
## Input

- PDBFILE:     PDBFILE
- HESSFILE:    Covariance (Inverse Hessian) Matrix in a [NxN] ascii format
- RESNUMS:     Chain + Residues number in the pdb, e.g. A15 B21

## Output Files

* Structure used for DFI: dfi-dficolor.pdb
* Master DFI: dfianalysis.csv

---

# Developers

- Avishek Kumar avishek@asu.edu

---

# References

1. Gerek, Z. N., Keskin, O. & Ozkan, S. B. Identification of specificity and promiscuity of PDZ domain interactions
through their dynamic behavior. Proteins 77, 796–811 (2009).
2. Gerek, Z. N. & Ozkan, S. B. Change in allosteric network affects binding affinities of PDZ domains: analysis
through perturbation response scanning. PLoS computational biology 7, e1002154 (2011).
3. Glembo, T. J., Farrell, D. W., Gerek, Z. N., Thorpe, M. & Ozkan, S. B. Collective dynamics differentiates
functional divergence in protein evolution. PLoS computational biology 8, e1002428 (2012).
4. Nevin Gerek, Z., Kumar, S. & Banu Ozkan, S. Structural dynamics flexibility informs function and evolution
at a proteome scale. Evolutionary applications 6, 423–33 (2013).
5. Woodrum, B. W., Maxwell, J. D., Bolia, A., Ozkan, S. B. & Ghirlanda, G. The antiviral lectin cyanovirin-N:
probing multivalency and glycan recognition through experimental and computational approaches. Biochemical Society transactions 41, 1170–6 (2013).
6. Bolia, A., Gerek, Z. N. & Ozkan, S. B. BP-Dock: a flexible docking scheme for exploring protein-ligand
interactions based on unbound structures. Journal of chemical information and modeling 54, 913–25 (2014).
7. Bolia, A. et al. A Flexible Docking Scheme Efficiently Captures the Energetics of Glycan-Cyanovirin
Binding. Biophysical Journal 106, 1142–1151 (2014).
8. Zou, T., Risso, V. A., Gavira, J. A., Sanchez-Ruiz, J. M. & Ozkan, S. B. Evolution of Conformational
Dynamics Determines the Conversion of a Promiscuous Generalist into a Specialist Enzyme. Molecular
biology and evolution (2014). doi:10.1093/molbev/msu281
9. Li, Z., Bolia, A., Ozkan, S. B., Ghirlanda, G. & Margulis, C. J. How hinge flexibility alters glycan
binding affinity of cyanovirin. in preparation (2014).
10. Bolia, A. & Ozkan, S. B. Adaptive BP-Dock: An Induced Fit Docking Approach for Full Receptor Flexibility.
Journal of Chemical Information and Modeling 56, 734–746 (2016).
11. Kim, H. et al. A Hinge Migration Mechanism Unlocks the Evolution of Green-to-Red Photoconversion in GFP-like
Proteins. Structure 23, 34–43 (2015).
12. Li, Z. et al. A Rigid Hinge Region Is Necessary for High-Affinity Binding of Dimannose to Cyanovirin and
Associated Constructs. Biochemistry 54, 6951–6960 (2015).
13. Butler, B. M., Gerek, Z. N., Kumar, S. & Ozkan, S. B. Conformational dynamics of nonsynonymous variants at protein
interfaces reveals disease association: The Role of Dynamics in Neutral and Damaging nsSNVs. Proteins: Structure,
Function, and Bioinformatics 83, 428–435 (2015).
14. Kumar, A., Butler, B. M., Kumar, S. & Ozkan, S. B. Integration of structural dynamics and molecular evolution
via protein interaction networks: a new era in genomic medicine. Current Opinion in Structural Biology 35, 135–142 (2015).
15. Kumar, A., Glembo, T. J. & Ozkan, S. B. The Role of Conformational Dynamics and Allostery in the Disease
Development of Human Ferritin. Biophysical Journal (2015). doi:10.1016/j.bpj.2015.06.060

---

# Funding and Acknowledgements

The DFI library was develped by Avishek Kumar (avishek@asu.edu) at Arizona State
University in the Ozkan Group under the following grants NIH-U54GM094599-05 and
NIH-R21LM011941-02. It is licensed under the BSD 3-clause license.

---
