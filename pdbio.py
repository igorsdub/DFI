#!/usr/bin/env python
#For reading and writing pdb files 

import numpy as np 
import sys

def read_pdb(filename):
    """
    DEPRICATED FUNCTION 
    function: 
        read_pdb, reads the alpha carbons of a model 
    arguments:
        pdbfile,string filename of a pdb file
    returns:
        Pos, a (N,3) dimensional NumPy array of positions
        ResNames, a N-length list containing the three-letter code for
           the amino acid type. 
    """
    with open(filename) as PDB:
        ResNames=[]
        Pos=[]
        for line in PDB:
            if line.startswith('TER') or line.startswith('ENDMDL'):
                break
            if line.startswith('ATOM') and line[13:16].strip(' ')=='CA':
                ResNames.append(line[17:20])
                Pos.append([line[31:38],line[39:46],line[47:54]])
    return np.array(Pos,dtype=float), ResNames 


def pdb_reader(filename,ATOMS,CAonly=False,noalc=True,chainA=False):
    readatoms=0
    with open(filename) as pdb:
        for line in pdb:
            if line.startswith('ENDMDL'):
                print "MULTIPLE MODELS...USING MODEL1"
                return 
            
            if line.startswith('ATOM'):
                record = line[:6]
                atom_index = line[7:11]
              
                atom_name = line[13:16]
                if CAonly and not(atom_name=='CA '):
                    continue 
                
                alc = line[16] #alternate location
                if noalc and not((alc==' ' or alc=='A')):
                    continue 

                res_name = line[17:20]
                
                chainID=line[21]
                if chainA and not(chainID=='A'):
                    continue 
                   
                res_index = line[23:26]
                insert_code = line[26]
                x = line[31:38]
                y = line[39:46]
                z = line[47:54]
                occupancy = line[55:60]
                temp_factor = line[61:66]
                ATOMS.append( ATOM(line[:6], line[7:11], line[13:16], line[16], line[17:20], line[21], line[23:26],
                                   line[26], line[31:38], line[39:46], line[47:54], line[55:60], line[61:66]) )
                readatoms+=1
    print "Read %d atoms from the %s"%(readatoms,filename)

def pdb_writer(ATOMS,msg="HEADER  frodaN unfolding target\n",filename="out.pdb",modelnum=1,atomoffset=0,residueoffset=0,mode="w"):
    
    with open(filename,mode) as pdb:
        pdb.write(msg)
        pdb.write("MODEL %d\n"%modelnum)
        pdb.write("PARENT N/A\n")
        for atom in ATOMS:
            record = atom.record
            atom_index=atom.atom_index + atomoffset 
            atom_name = atom.atom_name
            alc = atom.alc 
            res_name = atom.res_name
            chainID=atom.chainID
            res_index = atom.res_index + residueoffset 
            iCode = atom.insert_code
            x = atom.x
            y = atom.y
            z = atom.z
            occupancy = 1.00
            temp_factor = atom.temp_factor
            pdb.write("ATOM  %(atom_index)5d %(atom_name)4s%(alc)1s%(res_name)-3s %(chainID)1s%(res_index)4d%(iCode)1s   %(x)8.3f%(y)8.3f%(z)8.3f%(occupancy)6.2f%(temp_factor)6.2f\n" % vars())
        pdb.write("TER\n")
        pdb.write("END\n")
    print "Wrote out to file, %s"%filename
        
                
class ATOM:
    def __init__(self,record,atom_index,atom_name,alc,res_name,chainID,res_index,insert_code,x,y,z,occupancy,temp_factor):
        self.record=str(record)
        self.atom_index = int(atom_index)
        self.atom_name = str(atom_name)
        self.alc = alc
        self.res_name = str(res_name)
        self.chainID = str(chainID)
        self.res_index = int(res_index)
        self.insert_code = str(insert_code)
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.occupancy = float(occupancy) 
        self.temp_factor = float(temp_factor)


                

        
