#!/usr/bin/env python
"""
USAGE: ./PdbMunger PDBID 

Extracts the first chain of the PDB and Outputs to a file. 
TODO: Ouput only the alpha carbons of the pdbfile 
TODO: Ouput input parameters for the DFI Code. 

PDBID PDBID of the pdb file

"""

import sys

if len(sys.argv) < 2:
    print __doc__
    exit()

import os
from Bio import PDB

class PDBChainSlice:
    def __init__(self, out_dir=None):
        """ Create the writer and reader and destination for the fil. """
        self.parser = PDB.PDBParser()
        self.writer = PDB.PDBIO()
        if out_dir is None:
            out_dir=os.getcwd()
        self.out_dir = out_dir

    def make_pdb(self, pdb_path, chain_letters, overwrite=False, struct=None):
        """ Create a new PDB file containing only the specified chains.

        Returns the path to the created file.

        :param pdb_path: full path to the crystal structure
        :param chain_letters: iterable of chain characters (case insensitive)
        :param overwrite: write over the output file if it exists
        """
        chain_letters = [chain.upper() for chain in chain_letters]

        # Input/output files
        (pdb_dir, pdb_fn) = os.path.split(pdb_path)
        pdb_id = pdb_fn[3:7]
        out_name = "pdb%s_%s.pdb" % (pdb_id, "".join(chain_letters))
        out_path = os.path.join(self.out_dir, out_name)
        print "OUT PATH:",out_path
        plural = "s" if (len(chain_letters) > 1) else ""  # for printing

        # Skip PDB generation if the file already exists
        if (not overwrite) and (os.path.isfile(out_path)):
            print("Chain%s %s of '%s' already extracted to '%s'." %
                    (plural, ", ".join(chain_letters), pdb_id, out_name))
            return out_path

        print("Extracting chain%s %s from %s..." % (plural,
                ", ".join(chain_letters), pdb_fn))

        # Get structure, write new file with only given chains
        if struct is None:
            struct = self.parser.get_structure(pdb_id, pdb_path)
        self.writer.set_structure(struct)
        self.writer.save(out_path, select=SelectChains(chain_letters))

        return out_path

class SelectChains(PDB.Select):
    """ Only accept the specified chains when saving. """
    def __init__(self, chain_letters):
        self.chain_letters = chain_letters

    def accept_chain(self, chain):
        return (chain.get_id() in self.chain_letters)


def ExtractChain(PDBID,chainID='A'):
    """ Download the PDB and extract the chain (default A) """
    pdbList = PDB.PDBList()
    slicer = PDBChainSlice()


    pdb_fn = pdbList.retrieve_pdb_file(PDBID)
    
    return slicer.make_pdb(pdb_fn,chainID)

def extractA_CA(pdbcode,Verbose=False):
    """ Gets the pdb from the databank and then extracts chain A and the alpha carbons """
    
    outpath = ExtractChain(sys.argv[1])
    print "Outpath: %s"%outpath
    (pdb_dir, pdbfile) = os.path.split(outpath)
   
    pdbid=pdbfile.strip('_A.pdb')
    if Verbose:
        print pdbid
        
    os.system('grep "^ATOM.*CA" ' + pdbfile + ' > ' + pdbid+'_A_CA.pdb')
    pathtofile = os.path.join(os.getchwd(),"%s_A_CA.pdb"%(pdbid))
    return pathtofile  

if __name__ == "__main__":
    print extractA_CA(sys.argv[1])
