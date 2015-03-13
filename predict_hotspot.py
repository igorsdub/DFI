#!/usr/bin/env python 
print __doc__

import math,operator
import sys, os

def calcSASA(naccessToolPath, naccessOutputPath, pdb_file):
    pdb_ID = pdb_file
    pdb_file = "temp/%s.pdb" % (pdb_file)
    cmd = '%s/naccess %s' % (naccessToolPath, pdb_file)
    os.system(cmd)
    os.system("mv %s.rsa %s/%s.naccess" % (pdb_ID, naccessOutputPath, pdb_ID))
    os.system("rm *.log *.asa")


def parseChainAtoms(pdbPath, pdbID, chain):
    chain_coor = []
    # open the pdb file
    fileDir = "%s/%s.pdb" % (pdbPath,pdbID)
    inputFile = open(fileDir, "r")
    # open the file to write atoms
    flag = False
    line = inputFile.readline()
    outputFile = open("temp/%s%s.pdb" %(pdbID, chain),"w")
    while line:
        # record the atom if it belongs to the chain or this is a complex target
        if line.startswith("ATOM") and (line[21] == chain):
            outputFile.writelines(line)
            chain_coor.append(line)
            flag = True
        # if desired chain found or complex is recorded, break
        if flag == True and (line.startswith("ENDMDL") or (line.startswith("TER") and chain != " ")):
            break
        line = inputFile.readline()
    outputFile.writelines('END\n')
    outputFile.close()
    # end of while
    return chain_coor

def twoChainWrite(pdbPath,pdbID, ch1, ch2):
    chain1 = parseChainAtoms(pdbPath, pdbID, ch1)
    chain2 = parseChainAtoms(pdbPath, pdbID, ch2)
    complexfile = open('temp/%s%s%s.pdb' % (pdbID, ch1,ch2), 'w')
    for chline in chain1:
        complexfile.writelines(chline)
    complexfile.writelines('TER\n')
    for chline in chain2:
        complexfile.writelines(chline)
    complexfile.writelines('END\n')
    complexfile.close()

#    
#    main 
#
##################
def RunNaccess(pdbPath, interfaceID, naccessToolPath, naccessOutputPath):
    pdb_id = interfaceID[0:4]
    ch1 = interfaceID[4]
    ch2 = interfaceID[5]
    ch1_data =  parseChainAtoms(pdbPath, pdb_id, ch1) #PDBChain.getChain(pdbFile, ch1)
    if len(ch1_data) == 0:
        print "ERROR: Can not find chain", ch1, "in PDB file", pdb_id
        return
    ch2_data = parseChainAtoms(pdbPath, pdb_id, ch2)
    if len(ch2_data) == 0:
        print "ERROR: Can not find chain", ch2, "in PDB file", pdb_id
        return
    pf_data = twoChainWrite(pdbPath, pdb_id, ch1, ch2)
    rsa_file = ['']
    for structure in [pdb_id+ch1, pdb_id+ch2, pdb_id+ch1+ch2]:
        calcSASA(naccessToolPath, naccessOutputPath, structure)
    os.system("rm temp/*.asa temp/*.log")





def residueDict(pdbPath,pdbID):
    file = open(pdbPath + pdbID + ".pdb","r")
    resDict = {}
    while 1:
        line = file.readline()
        if line == "":
            break
        if line.startswith("ENDMDL"):
            break
        if line.startswith("ATOM"):
            resno = line[22:26].replace(" ","")
            chain = line[21]
            resid = line[17:20]
            resDict[resno+chain] = resid
    file.close()
    return resDict


def chainFetch(interfaceID):
    pdbID = interfaceID[0:4]
    chain1 = interfaceID[4]
    chain2 = interfaceID[5]
    chains = set()
    file = open(pdbPath + pdbID+'.pdb', 'r')
    while 1:
        line = file.readline()
        if line == '':
            break
        if line.startswith('ENDMDL'):
            break
        if line.startswith('ATOM'):
            chains.add(line[21])
    file.close()
    chains=list(chains)
    
    if chains.count(chain1) ==0 and chains.count(chain2) ==1:
	print "Chain %s does not exist in the pdb file! \n" %(chain1)
	return -1
    elif chains.count(chain1) ==1 and chains.count(chain2) ==0:
	print "Chain %s does not exist in the pdb file! \n" %(chain2)
	return -1
    elif chains.count(chain1) ==0 and chains.count(chain2) ==0:
	print "Chain %s and Chain %s do not exist in the pdb file! \n" %(chain1,chain2)
	return -1
    
    return chains

def three2one(residue):
    if residue ==   "ALA": return "A"
    elif residue == "CYS": return "C"
    elif residue == "ASP": return "D"
    elif residue == "GLU": return "E"
    elif residue == "PHE": return "F"
    elif residue == "GLY": return "G"
    elif residue == "HIS": return "H"
    elif residue == "ILE": return "I"
    elif residue == "LYS": return "K"
    elif residue == "LEU": return "L"
    elif residue == "MET": return "M"
    elif residue == "ASN": return "N"
    elif residue == "PRO": return "P"
    elif residue == "GLN": return "Q"
    elif residue == "ARG": return "R"
    elif residue == "SER": return "S"
    elif residue == "THR": return "T"
    elif residue == "VAL": return "V"
    elif residue == "TRP": return "W"
    elif residue == "TYR": return "Y"
    else: return "X"

## Distance calculation method
def distanceCalculation(coordinates1, coordinates2):
    x1 = coordinates1[0]
    y1 = coordinates1[1]
    z1 = coordinates1[2]
    x2 = coordinates2[0]
    y2 = coordinates2[1]
    z2 = coordinates2[2]
    distance = math.sqrt(((x2-x1)**2)+((y2-y1)**2)+((z2-z1)**2))
    return distance

## Contacting residues extraction
def CoordinatesFetch(pdbPath, pdbID, chain1, chain2):
##    file = open("PairPotentialsDict","r")
##    pairPotentialDict = pickle.load(file)
##    file.close()
    chain1Atoms = []
    chain2Atoms = []
##    d = {'H':1.20, 'C':1.70, 'N':1.55, 'O':1.52, 'F':1.35, 'P':1.90, 'S':1.85}
    VDW_RADII = { 'C':1.76, 'N':1.65, 'O':1.40, 'CA':1.87, 'H':1.20, 'S':1.85, 'CB':1.87, 'CZ':1.76, 'NZ':1.50, 'CD':1.81, 'CE':1.81, 'CG':1.81, 'C1':1.80, 'P':1.90 }

    VDW_RADII_EXTENDED = { 'ALA': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'OXT': 1.40  },\
                           'ARG': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.87, 'CD': 1.87, 'NE': 1.65, 'CZ': 1.76, 'NH1': 1.65, 'NH2': 1.65, 'OXT': 1.40 },\
                           'ASP': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.76, 'OD1': 1.40, 'OD2': 1.40, 'OXT': 1.40  },\
                           'ASN': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.76, 'OD1': 1.40, 'ND2': 1.65, 'OXT': 1.40 },\
                           'CYS': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'SG': 1.85, 'OXT': 1.40 },\
                           'GLU': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.87, 'CD': 1.76, 'OE1': 1.40, 'OE2': 1.40 , 'OXT': 1.40 },\
                           'GLN': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.87, 'CD': 1.76, 'OE1': 1.40, 'NE2': 1.65, 'OXT': 1.40  },\
                           'GLY': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'OXT': 1.40 },\
                           'HIS': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.76, 'ND1': 1.65, 'CD2': 1.76, 'CE1': 1.76, 'NE2': 1.65, 'OXT': 1.40  },\
                           'ILE': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG1': 1.87, 'CG2': 1.87, 'CD1': 1.87 , 'OXT': 1.40 },\
                           'LEU': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.87, 'CD1': 1.87, 'CD2': 1.87, 'OXT': 1.40 },\
                           'LYS': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.87, 'CD': 1.87, 'CE': 1.87, 'NZ': 1.50 , 'OXT': 1.40 },\
                           'MET': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.87, 'SD': 1.85, 'CE': 1.87, 'OXT': 1.40 },\
                           'PHE': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.76, 'CD1': 1.76, 'CD2': 1.76, 'CE1': 1.76, 'CE2': 1.76, 'CZ': 1.76, 'OXT': 1.40  },\
                           'PRO': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.87, 'CD': 1.87, 'OXT': 1.40  },\
                           'SER': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'OG': 1.40, 'OXT': 1.40 },\
                           'THR': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'OG1': 1.40, 'CG2': 1.87, 'OXT': 1.40 },\
                           'TRP': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.76, 'CD1': 1.76, 'CD2': 1.76, 'NE1': 1.65, 'CE2': 1.76, 'CE3': 1.76, 'CZ2': 1.76, 'CZ3': 1.76, 'CH2': 1.76, 'OXT': 1.40  },\
                           'TYR': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.76, 'CD1': 1.76, 'CD2': 1.76, 'CE1': 1.76, 'CE2': 1.76, 'CZ': 1.76, 'OH': 1.40, 'OXT': 1.40 },\
                           'VAL': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG1': 1.87, 'CG2': 1.87, 'OXT': 1.40  }, 'ASX': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.76, 'AD1': 1.50, 'AD2': 1.50 }, 'GLX': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.76, 'CD': 1.87, 'AE1': 1.50, 'AE2': 1.50, 'OXT': 1.40  },\
                           'ACE': { 'C': 1.76, 'O': 1.40, 'CA': 1.87 }, 'PCA':VDW_RADII, 'UNK':VDW_RADII }

   
    file = open(pdbPath + pdbID + ".pdb","r")
    atomList = []
    contact = set()
    left = set()
    right = set()
    while 1:
        line = file.readline()
        if line == "":
            break
        if line.startswith("ENDMDL") or line.startswith("ENDMODEL"):
            break
        if line.startswith("ATOM"):
            if line[26]==" ":
                atom = line[13:16].replace(" ","")
                resid = line[17:20]
                resno = line[22:26].replace(" ","")
##                print resno
                chain = line[21]
                atomType = atom[0]
                coordinates = [float(line[27:38]), float(line[38:46]), float(line[46:54])]

                if chain == chain1:
                    chain1Atoms.append([resno, coordinates, atomType, resid, atom])
                if chain == chain2:
                    chain2Atoms.append([resno, coordinates, atomType, resid, atom])

##three2one(resid)+resno+chain
    for atom1 in chain1Atoms:
        for atom2 in chain2Atoms:
            #print atom1[3], atom2[3]
            d1 = VDW_RADII_EXTENDED[atom1[3]]
            d2 = VDW_RADII_EXTENDED[atom2[3]]
            ## print atom1, atom2
##            if atom1[4] == "OXT":
##                print atom1[3], atom1[4]
##
##            if atom2[4] == "OXT":
##                print atom2[3], atom2[4]
            if three2one(atom1[3])!= 'X' and atom1[4].startswith('H')== False and three2one(atom2[3])!= 'X'and atom2[4].startswith('H')== False:
                #print atom1[4], atom2[4]
                try:
                    cutoff = d1[atom1[4]]+d2[atom2[4]]+ 0.5
                    #print atom1, atom2
                    distance = distanceCalculation(atom1[1], atom2[1])
                    #print cutoff
                    if distance <= cutoff:
                        #print distance
                        cont = [atom1[0]+chain1,atom2[0]+chain2]
                        cont.sort()
                        left.add(three2one(atom1[3])+atom1[0]+chain1)
                        right.add(three2one(atom2[3])+ atom2[0]+chain2)
                        residList = [three2one(atom1[3]),three2one(atom2[3])]
                        residList.sort()
                except KeyError:
                    continue
                    #print atom1, atom2
                    continue
##                invDist = pairPotentialDict[residList[0]+"-"+residList[1]]
##                contact.add(cont[0]+'\tpp\t'+cont[1]+ '\t%5.2f' % invDist)
    file.close()
    left = list(left)
    left.sort()
    right = list(right)
    right.sort()
    if len(left) <= 0 and len(right)<=0 :
	print "There is no interface between Chain %s and Chain %s !\n" %(chain1, chain2)
	return -1
    else:
        return left, right
    
    
    ## Nearby residue extraction.
def nearbyResidues(pdbPath, pdbID, chainID, residues):
    nearbyResidue = []
    #~ file = open("PairPotentialsDict","r")
    #~ pairPotentialDict = pickle.load(file)
    #~ file.close()
    file = open(pdbPath + pdbID + ".pdb","r")
    CAset = {}
    while 1:
        line = file.readline()
        if line == "":
            break
        if line.startswith("ENDMDL"):
            break
        if line.startswith("ATOM"):
            atom = line[13:16].replace(" ","")
            resno = line[22:26].replace(" ","")
            chain = line[21]
            atomType = atom[0]
            resid = three2one(line[17:20])
            coordinates = [float(line[27:38]), float(line[38:46]), float(line[46:54])]
            if atom == 'CA':
                if chain == chainID:
                    CAset[resid+resno+chain] =  coordinates
    file.close()
    for residue in residues:
        for CA in CAset.keys():
            if CA != residue:
             
                resid1 = int(CA[1:-1])
                resid2 = int(residue[1:-1])
                if abs(resid1-resid2) > 1:
                    distance = distanceCalculation(CAset[residue], CAset[CA])
                    if distance <= 6.0:
                        #residList = [three2one(residue.split(".")[1]), three2one(CA.split(".")[1])]
                        #residList.sort()
                        #~ invDist = pairPotentialDict[residList[0]+"-"+residList[1]]
                        nearbyResidue.append(residue + '\tpp\t' + CA )
                if abs(resid1-resid2) == 1:
                    distance = distanceCalculation(CAset[residue], CAset[CA])
                    if distance <= 6.0:
                        #residList = [three2one(residue.split(".")[1]), three2one(CA.split(".")[1])]
                        #residList.sort()
                        #~ invDist = pairPotentialDict[residList[0]+"-"+residList[1]]  # * 10
                        nearbyResidue.append(residue + '\tpp\t' + CA )
                   # print residue.split(".")[0] + '\tpp\t' + CA.split(".")[0] + '\t%5.2f' % invDist
    return nearbyResidue
    
    
    ## PDB file generator
def PDBfileGenerator(interfaceID):
    pdbID = interfaceID[0:4]
    chain1 = interfaceID[4]
    chain2 = interfaceID[5]
    left, right = CoordinatesFetch(pdbPath, pdbID, chain1, chain2)
##    for item in contact:
##        print item
    chain1Residues = set()
    chain2Residues = set()
    for item in left:
        chain1Residues.add(int(item[1:-1]))
    for item in nearbyResidues(pdbPath, pdbID, chain1, left):
#	print "Item: "+item
        chain1Residues.add(int(item.split()[2][1:-1]))
    for item in right:
        chain2Residues.add(int(item[1:-1]))
    for item in nearbyResidues(pdbPath, pdbID, chain2, right):
        chain2Residues.add(int(item.split()[2][1:-1]))
    chain1Residues = list(chain1Residues)
    chain1Residues.sort()
    chain2Residues = list(chain2Residues)
    chain2Residues.sort()
##    print chain1Residues
##    print chain2Residues
    pdblines = set()
    newfile = open(pdbPath + pdbID + ".pdb","r")
    file = open(resultPath  + interfaceID + '.pdb', 'w')
    #~ print "interfacePdb opened"
    while 1:
        line = newfile.readline()
        if line == "":
            break
        if line.startswith("ENDMDL"):
            break
        if line.startswith("ATOM"):
            resno = line[22:26].replace(" ","")
            chain = line[21]
            if chain == chain1:
                for item in chain1Residues:
                    if resno == str(item):
                        file.writelines(line)
            if chain == chain2:        
                for item in chain2Residues:
                    if resno == str(item):
                        file.writelines(line)
                    
    newfile.close()
    file.close()
    return 1

######################    PAIR POTENTIAL    ################################
############################################################################

###################################### Center of Mass #########################################################################
## Specify the residue atoms of interest to extract the center of mass coordinates. Dump these coordinates into a dictionary,
## where keys are residue ID (residue name + residue number + chain), the values are x,y,z coordinates.
## Use this function in "contactingResidues" function to extract contacting pairs.
    
def center_of_mass(pdbName, chain1, chain2):
    alphaCarbons = Ca_CoordinatesFetch_BothChain(pdbName, chain1, chain2)
    ## dictionary of interested atoms which are sent by okeskin.
    dictAtom = {'ACE':['CA'], 'ALA':['CB'], 'GLY':['CA'], 'SER':['OG'], 'ASN':['OD1', 'ND2'], 'ASP':['OD1', 'OD2'], 'CYS':['SG'], 'PRO':['CB', 'CG', 'CD'], \
           'THR':['OG1'], 'VAL':['CG1', 'CG2'], 'ARG':['NE', 'NH1', 'NH2'], 'GLN':['OE1', 'NE2'], 'GLU':['OE1', 'OE2'], 'HIS':['CG', 'ND1', 'CD2', 'CE1', 'NE2'], 'ILE':['CD1'], \
           'LEU':['CD1', 'CD2'], 'LYS':['NZ'], 'MET':['SD'], 'PHE':['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'], 'TRP':['CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'], 'TYR':['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH']}
    centerCoordinates = {}
    for residues in alphaCarbons.keys():
        resname = residues[0:3]
        allAtoms = alphaCarbons[residues]
        x_total = 0
        y_total = 0
        z_total = 0
        coordinateDict = {}
        for ind_atom in allAtoms:
            coordinateDict[ind_atom[0]] = ind_atom[1]
        atoms = dictAtom[resname]
        try:
            for item in atoms:
                
                coordinates = coordinateDict[item]
                x = coordinates[0]
                y = coordinates[1]
                z = coordinates[2]
                x_total += x
                y_total += y
                z_total += z
    ##        print residues, residues[3:], residues[0:]
            real_coordinates = [x_total/len(atoms), y_total/len(atoms), z_total/len(atoms)]
            centerCoordinates[pdbName+three2one(resname)+residues[3:]] = real_coordinates
        except KeyError:
            # print pdbName, item, residues[3:]
            continue
    return centerCoordinates
###############################################################################################################################



####################################### Atom Fetch ############################################################################
## This function extracts all atoms and their coordinates for given pdb name and chain names.
## It returns the dictionary of atoms for individual residues in the protein chains.
## This function is used by "center_of_mass"

def Ca_CoordinatesFetch_BothChain(pdbName, chain1, chain2):
    alphaCarbons = {}
    file = open(pdbPath + pdbName + ".pdb","r")
    atomList = []
    while 1:
        line = file.readline()
        if line == "":
            break
        if line.startswith("ENDMDL"):
            break

        if line.startswith("ATOM"):
            atom = line[13:16].replace(" ","")
            resid = line[17:20]
            resno = line[22:26].replace(" ","")
            chain = line[21]
            coordinates = [float(line[27:38]), float(line[38:46]), float(line[46:54])]
            if chain == chain1 or chain == chain2:
                if line[26] == ' ':
                    try:
##                        print 
                        temp_list = alphaCarbons[resid+resno+chain]
##                        print  , resid+resno+chain
                        temp_list.append([atom, coordinates])
                        alphaCarbons[resid+resno+chain] = temp_list
                    except KeyError:
                        alphaCarbons[resid+resno+chain] = [[atom, coordinates]]
    file.close()
    return alphaCarbons
##############################################################################################################################

  
import math


######################################### Contacting residues according to center of mass ###################################
## This function extracts contacting residues each from one chain. The cutoff is defined as 7 Angstrom.
## Contacting pairs of each residue is dumped in a dictionary of each individual residue.
## This dictionary will be used to calculate delta pair potentials for each residue.

def contactingResidues(pdbName, chain1, chain2):
    centeredResidues = center_of_mass(pdbName, chain1, chain2)
    contactDict = {}
    for res1 in centeredResidues.keys():
        for res2 in centeredResidues.keys():
            coor1 = centeredResidues[res1]
            coor2 = centeredResidues[res2]
            if res1 != res2:
                if res1[-1] != res2[-1]:
                    distance = distanceCalculation(coor1, coor2)
                    if distance <= 7.0:
                        try:
                            temp = contactDict[res1]
                            temp.append(res2)
                            contactDict[res1] = temp
                        except KeyError:
                            contactDict[res1] = [res2]
                if res1[-1] == res2[-1]:
                    if abs(int(res1[5:len(res1)-1])-int(res2[5:len(res2)-1])) > 3:
                        distance = distanceCalculation(coor1, coor2)
                        if distance <= 7.0:
                            try:
                                temp = contactDict[res1]
                                temp.append(res2)
                                contactDict[res1] = temp
                            except KeyError:
                                contactDict[res1] = [res2]

    return contactDict
#############################################################################################################################


def contactPotentials(pdbName, chain1, chain2):
    contactDict = contactingResidues(pdbName, chain1, chain2)
    ddPP = {}
    allResidues = Ca_CoordinatesFetch_BothChain(pdbName, chain1, chain2)
    matrixPotentials = PairPot()
    for residues in contactDict.keys():
        totalContact = 0
        res1 = residues
        for res2 in contactDict[res1]:
            aa1 = res1[4]
            aa2 = res2[4]
            if aa1 != 'X' and aa2!= 'X':
                temp = [aa1, aa2]
                temp.sort()
                contact = temp[0] + "-" + temp[1]
                totalContact += matrixPotentials[contact]
        if len(contactDict[res1])!= 0:    
            ddPP[res1] = totalContact
        if len(contactDict[res1])== 0:
            ddPP[res1] = 0.0
    for item in allResidues.keys():
        try:
            ddPP[pdbName+three2one(item[0:3])+item[3:]]
        except KeyError:
            ddPP[pdbName+three2one(item[0:3])+item[3:]] = 0.0
    return ddPP

def PairPot():
    file = open('true_solvent_mediated_pairpotentials.list','r')
##    file = open('true_solvent_mediated_pairpotentials.list','r')
    index = 0
    resDict = {}
    matrix = []
    while 1:
        line = file.readline()
        if line == "":
            break
        temp = line.strip().split()
        resDict[index] = temp[0]
        matrixEntry = temp[1:]
        matrix.append([])
        for i in range(len(matrixEntry)):
            matrix[index].append(float(matrixEntry[i]))

        index += 1
    PairPotentialDict = {}
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            if i<=j:#i<=j:(solvent mediated)#if i >= j:(residue mediated)
                tempRes = [resDict[i], resDict[j]]
                tempRes.sort()
                pair = tempRes[0] + '-' + tempRes[1]
                PairPotentialDict[pair] = matrix[i][j]
    return PairPotentialDict


######################   INDIVIDUAL ASA FETCH ##############################
############################################################################

def NaccessFileParserComplex(interface):
    ASAcomplex = {}
    file = open(naccessOutputPath +interface+".naccess","r")
    while 1:
        line = file.readline()
        if line == "":
            break
        if line.startswith("RES"):
            if line[13] == " ":
                #temp = line.strip().split()
                #resName = three2one(temp[1])
                #resPosition = temp[3]
                #chain = temp[2]
                #ASA = float(temp[4])
                resName = three2one(line[4:7])
                resPosition = line[9:13].replace(" ","")
                chain = line[8]
                temp = line[14:].strip().split()
                ASA = float(temp[0])
                ASAcomplex[resName+resPosition+chain] = ASA
    file.close()
    return ASAcomplex

def NaccessFileParserMonomer(interface):
    ASAmonomer = {}
    monomer1 = interface[0:5]
    monomer2 = interface[0:4]+interface[5]
    file = open(naccessOutputPath + monomer1+".naccess" ,"r")
    while 1:
        line = file.readline()
        if line == "":
            break
        if line.startswith("RES"):
            #temp = line.strip().split()
            #resName = three2one(temp[1])
            #resPosition = temp[3]
            #chain = temp[2]
            #ASA = float(temp[4])
            resName = three2one(line[4:7])
            resPosition = line[9:13].replace(" ","")
            chain = line[8]
            temp = line[14:].strip().split()
            ASA = float(temp[0])
            ASAmonomer[resName+resPosition+chain] = ASA
    file.close()

    file = open(naccessOutputPath +monomer2+".naccess","r")
    while 1:
        line = file.readline()
        if line == "":
            break
        if line.startswith("RES"):
            #temp = line.strip().split()
            #resName = three2one(temp[1])
            #resPosition = temp[3]
            #chain = temp[2]
            #ASA = float(temp[4])
            resName = three2one(line[4:7])
            resPosition = line[9:13].replace(" ","")
            chain = line[8]
            temp = line[14:].strip().split()
            ASA = float(temp[0])
            ASAmonomer[resName+resPosition+chain] = ASA
            #ASAmonomer[resName+resPosition+chain] = ASA
    file.close()
    return ASAmonomer

def StandardASADict():
    dStandard = {}
    file = open("standard.data","r")
    lines = file.readlines()
    for line in lines:
        temp = line.strip().split()
        dStandard[three2one(temp[3])] = float(temp[4])
    file.close()
    return dStandard

############################################################################
############################################################################


dStandard = StandardASADict()
def HSPredictionDriver(interface):
    file = open("%s%s.result" % (resultPath, interface), "w")
    newfile = open("%s%s.output" % (resultPath, interface), "w")
    jmolfile = open("%s%s.spt" % (resultPath, interface), "w")
    
    ASAmonomer = NaccessFileParserMonomer(interface)
    ASAcomplex = NaccessFileParserComplex(interface)
    pdbID = interface[0:4]
    chain1 = interface[4]
    chain2 = interface[5]
    PairPotential = contactPotentials(pdbID, chain1, chain2)
    contact = contactingResidues(pdbID, chain1, chain2)
    ##print ASAmonomer
    file.writelines("Residue Number_Residue Name_Chain_RelCompASA_RelMonomerASA_Potential_Prediction\n")
    newfile.writelines("Residue Number,Residue Name,Chain,RelCompASA,RelMonomerASA,Potential,Prediction\n")
    #newfile.writelines("------------------------------------------------------------------------------------------------------\n")

    [left, right] = CoordinatesFetch(pdbPath, pdbID, chain1, chain2)
    temp1 =[]
    temp2 =[]

    for item in left:
        #print item
        try:
            compASA = ASAcomplex[item]
	    monomerASA = ASAmonomer[item]
            potential = PairPotential[pdbID+item]
            RelCompASA = compASA*100.0/dStandard[item[0]]
	    RelMonomerASA = monomerASA*100.0/dStandard[item[0]]
    ##    print contact["1a0o"+item]
            if RelCompASA <= 20.0 and abs(potential) >= 18.0:
	        temp1.append([ int(item[1:len(item)-1]) , item[0] , item[-1] , RelCompASA , RelMonomerASA, abs(potential), "H"])
            else:
	        temp1.append([ int(item[1:len(item)-1]) , item[0] , item[-1] , RelCompASA , RelMonomerASA, abs(potential), "NH"])
        except KeyError:
            continue
    for item in right:
        #print item
        try:
            compASA = ASAcomplex[item]
            monomerASA = ASAmonomer[item]
            potential = PairPotential[pdbID+item]
            RelCompASA = compASA*100.0/dStandard[item[0]]
            RelMonomerASA = monomerASA*100.0/dStandard[item[0]]
    ##    print contact["1a0o"+item]
            if RelCompASA <= 20.0 and abs(potential) >= 18.0:
	        temp2.append([ int(item[1:len(item)-1]) , item[0] , item[-1] , RelCompASA , RelMonomerASA, abs(potential), "H"])
            else:
	        temp2.append([ int(item[1:len(item)-1]) , item[0] , item[-1] , RelCompASA ,RelMonomerASA, abs(potential), "NH"])
        except KeyError:
            continue
    #file.close()
    temp1 = sorted(temp1 , key=operator.itemgetter(0))
    temp2 = sorted(temp2 , key=operator.itemgetter(0))
    
    #sessionID=resultPath.split('/')[-3]
    
    jmolfile.writelines("zap;\n")
    jmolfile.writelines("load %s.pdb;\n" % pdbID)
    jmolfile.writelines("background white;\n")    
    jmolfile.writelines("wireframe off;\n")    
    jmolfile.writelines("restrict none;\n")    
    jmolfile.writelines("select *:%s;\n" %(chain1))       
    jmolfile.writelines("select selected or *:%s;\n" %(chain2))    
    jmolfile.writelines("trace on;\n")     
    jmolfile.writelines("select *:%s;\n" %(chain1))    
    jmolfile.writelines("color cyan;\n")  
    jmolfile.writelines("select *:%s;\n" %(chain2)) 
    jmolfile.writelines("color orange;\n")    
    jmolfile.writelines("select none;\n")   

    for item in temp1:
	if item[6]=="H":
	    jmolfile.writelines("select selected or %s:%s;\n" %(item[0],item[2]))  

    jmolfile.writelines("color red;\n")  
    jmolfile.writelines("select none;\n")   
    
    for item in temp2:
	if item[6]=="H":
	    jmolfile.writelines("select selected or %s:%s;\n" %(item[0],item[2]))  
	    
    jmolfile.writelines("color purple;\n")  
    jmolfile.writelines("select none;\n")   
    
    
    
    for item in temp1:
	file.writelines("%d_%s_%s_%6.2f_%6.2f_%6.2f_%s\n" %(item[0] , item[1] , item[2] , item[3] , item[4] , item[5] , item[6]))
	newfile.writelines("%d,%s,%s,%6.2f,%6.2f,%6.2f,%s\n" %(item[0] , item[1] , item[2] , item[3] , item[4] , item[5] , item[6]))
	jmolfile.writelines("select selected or %s:%s;\n" %(item[0],item[2]))  

	
    for item in temp2:
	file.writelines("%d_%s_%s_%6.2f_%6.2f_%6.2f_%s\n" %(item[0] , item[1] , item[2] , item[3] , item[4] , item[5] , item[6]))
	newfile.writelines("%d,%s,%s,%6.2f,%6.2f,%6.2f,%s\n" %(item[0] , item[1] , item[2] , item[3] , item[4] , item[5] , item[6]))
	jmolfile.writelines("select selected or %s:%s;\n" %(item[0],item[2]))  

    jmolfile.writelines("wireframe off;\n")    
    jmolfile.writelines("spacefill on;\n")  
    jmolfile.writelines("define bothChains selected;\n")    
    jmolfile.writelines("select bothChains and *:%s;\n" %(chain1))    
    jmolfile.writelines("define chainOne selected;\n")   
    jmolfile.writelines("select bothChains and *:%s;\n" %(chain2))   
    jmolfile.writelines("define chainTwo selected;\n")       
    
    jmolfile.close()
    file.close()
    newfile.close()
	    
    return 1

if len(sys.argv) == 6:
    interfaceID = sys.argv[1]
    pdbPath = sys.argv[2]
    naccessOutputPath = sys.argv[3]
    resultPath = sys.argv[4]
    naccessToolPath = sys.argv[5]
    ret = chainFetch(interfaceID)
    if  ret==-1:
	print "Try another structure\n"
    else:
        RunNaccess(pdbPath, interfaceID, naccessToolPath,naccessOutputPath)
        PDBfileGenerator(interfaceID)
        HSPredictionDriver(interfaceID)
else:
    print "usage: python predict_hotspot.py <interface> <pdbPath> <naccessOutputPath> <resultPath> <naccessToolPath>"

