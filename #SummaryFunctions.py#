import sys
import os
import shutil
import numpy
from numpy import linalg
from math import sqrt
import math
import urllib, urllib2
import sequence
import decimal
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio import Entrez
Entrez.email = "zgerek@asu.edu"
from Bio.Blast.Applications import NcbiblastpCommandline
from StringIO import StringIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from collections import deque
from decimal import *

workDir = os.curdir

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.PDBList import PDBList
from Bio.PDB.PDBIO import PDBIO 
from Bio.PDB.PDBIO import Select
from Bio.PDB.DSSP import DSSP

# this path should be changed!
DSSPPath='/Users/kumarlab/scripts/dssp/dsspcmbi'

class chain_select(Select):
	def accept_chain(self, chain):
		print(chain.get_id())
		if chain.get_id() == 'A':
			return 1
		else:
			return 0

class Chdir:         
      def __init__( self, newPath ):  
        self.savedPath = os.getcwd()
        os.chdir(newPath)

      def __del__( self ):
        os.chdir( self.savedPath )

 	
#class NotDisordered(Select):
#    def accept_atom(self, atom):
#        return not atom.is_disordered() or
#               atom.get_altloc()=='A'


def unique(items):
    found = set([])
    keep = []

    for item in items:
        if item not in found:
            found.add(item)
            keep.append(item)

    return keep

class FastaSeq:
	def __init__(self, name, sequence):
    		self.name = name
    		self.sequence = sequence

	def get_seqs(file):
    		items = []
    		index = 0
    		for line in file:
        		if line.startswith(">"):
            			if index >= 1:
                			items.append(aninstance)
            			index+=1
            			name = line[:-1]
            			seq = ''
            			aninstance = FastaSeq(name, seq)
        		else:
            			seq += line[:-1]
            			aninstance = FastaSeq(name, seq)
    		items.append(aninstance)

   		return items

def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


##retrieves pdb file##
def getPDBFiles(pdb_id):
    pdb_file_name = pdb_id + ".pdb"
    pdb_file = open(pdb_file_name, "w")
    url = 'http://www.rcsb.org/pdb/files/%s.pdb' % pdb_id
    pdb = urllib.urlopen(url).read()
    pdb_file.write(pdb)	
    pdb_file.close()
    
##--------------------------------##
##creates a file for CA##
##--------------------------------##
def createCAfile(pdb_id):
        cwd = os.getcwd()
        print cwd
        monomer_pdbFile =cwd+'/'+pdb_id+'.pdb'
        j = 0
#        readpdb_file = file(monomer_pdbFile,"rU").read()
#        writepdb_file = open(pdb_id + '-CA.pdb', "w")
#        pdblines = readpdb_file.split('\n')
#        OutPdb=deque([])
#        for line in pdblines:
#            if line.startswith("ATOM"):
#                if line[12:16] == " CA ":
#                      OutPdb.append(line)
#                      i+= i+1
#                      
#        totalCANumber = len(OutPdb)
#
#        while (i != j):
#            writepdb_file.write(OutPdb.popleft())
#            writepdb_file.write('\n')
#            j+= j+1
#        writepdb_file.write('END\n')
#        writepdb_file.close()
        
        p=PDBParser()

        TWO = Decimal(10) ** -2
        Decimal(3)
        structure = p.get_structure('InPdb', monomer_pdbFile)
        readpdb_file = open(monomer_pdbFile, "rb")
        writepdb_file = open(pdb_id + '-CA.pdb', "w")          
        for model in structure:
            for chain in model: 
               for residue in chain:            
                        for atom in residue:
                            if(atom==residue['CA']):
                                s = '{0:6}{1:5}{2:5}{3:4}{4:2}{5:8}{6:8}{7:8}{8:8}{9:6}{10:6}\n'.format('ATOM', atom.serial_number, atom.name.center(5), residue.resname.rjust(4),chain.get_id().rjust(2),str(residue.id[1]).rjust(4), Decimal(str(atom.coord[0]).rjust(9)), Decimal(str(atom.coord[1])), Decimal(str(atom.coord[2])), Decimal(str(atom.occupancy)).quantize(TWO), Decimal(str(atom.bfactor)).quantize(TWO))
                                j += 1
                                writepdb_file.write(s)
        writepdb_file.write('END\n')
        writepdb_file.close()

        totalCANumber = j 
        return totalCANumber

##---------------------------------------------##
##Replaces b-factor with %dfi
##---------------------------------------------##
def replaceBValues(File, pdbid, s2perc):
        length=len(s2perc)
        print length
        nMax = 5000;
      	ResNum = numpy.zeros(nMax,float)
        i = 0
        icount = 0
        readpdb_file = open(File, "rb")
	writedfi_file = open(pdbid+'-dfi.pdb', "w")
	

        for line in readpdb_file:
                if line.startswith("ATOM"):
                        ThisResNum = int(line[22:29])
                        if line[12:16] == " CA ":
                              ResNum[icount] = int(line[22:29])
                              icount +=1
        totCA = icount
        
      	readpdb_file.close()
        readpdb_file = open(File, "rb")
	for line in readpdb_file:
              if line.startswith("ATOM"):
                  ThisResNum = int(line[22:29])
                  for i in range(0,totCA+1):
                      if ThisResNum == ResNum[i]:
                                temp = float(s2perc[i])
                                line = line[:60] + "  " + str(decimal.Decimal('%.3f' % temp)) + line[66:]
                                writedfi_file.write(line)
              
	writedfi_file.write('END\n')
	writedfi_file.close()


#######################
### Blast with NPid ###
#######################
	

def npBLAST(NP_id):
	out_file= open(NP_id + ".fasta", "w")
	net_handle = Entrez.efetch(db="nucleotide", id=NP_id, rettype="fasta")
	out_file.write(net_handle.read())
	out_file.close()
	result_handle = NCBIWWW.qblast("blastp", "pdb", NP_id) 
	save_file = open(NP_id + "_blast.xml", "w")
	save_file.write(result_handle.read())
	save_file.close()
	result_handle.close()

#npBLAST("NP_000289")	

def seqBLAST(NP_id):
	handle = open(NP_id + ".fasta", "rU")
	sequenceHandle = SeqIO.parse(handle, "fasta")
	for record in sequenceHandle:
		sequence = str(record.seq)
	handle.close()
	result_handle = NCBIWWW.qblast("blastp", "pdb", sequence) 
	save_file = open(NP_id + "_blast.xml", "w")
	save_file.write(result_handle.read())
	save_file.close()
	result_handle.close()
	
#seqBLAST("NP_000289")

def parseBlastFile(NP_id):
	result_handle = open(NP_id + "_blast.xml")        
	blast_record = NCBIXML.read(result_handle)
	E_VALUE_THRESH = 1E-25
	out_file= open(NP_id+'.out',"w")
	print >> out_file,'#'+NP_id
	sequencequeryLength = blast_record.query_length
	for alignment in blast_record.alignments:
		for hsp in alignment.hsps:
			if hsp.expect < E_VALUE_THRESH:
				first = float(hsp.identities)
				second = len(hsp.query)
				identity = 100*float(first/second)
				identity = int(round(identity,0))
				coverage = 100*float(first/sequencequeryLength)
				coverage = int(round(coverage,0))
				line1=alignment.title
				b=line1.split('|')
				out_file.write(str(b[3]) +"\t")
#				print >> out_file,b[3],'\t'
				line2=b[4]
				c=line2.split()
#				print >> out_file,c[0],'\t'
				out_file.write(str(c[0]) +"\t")
				out_file.write(str(coverage) +"\t")
				out_file.write(str(identity) +"\t")
				out_file.write(str(hsp.expect) +"\t")
				out_file.write(str(hsp.query_start) +"\t")
				out_file.write(str(hsp.query_end) +"\t")
				
				out_file.write(str(hsp.sbjct_start) +"\t")
				out_file.write(str(hsp.sbjct_end) +"\n")




					
#parseBlastFile("NP_000289")

############################
### Blast with Uniprotid ###
############################				
				
def unBLAST(uniprot_id):

        #retrieves gene sequence
        url = 'http://uniprot.org/uniprot/'+ uniprot_id + '.fasta' 
        
        request = urllib2.Request(url)
        response = urllib2.urlopen(request)
        page = response.read(200000)
        out_file= open(uniprot_id + ".fasta", "w")
	out_file.write(page)
	out_file.close()
	

#------ finding matching pdb

def matchPDB(NP_id, minQueryCoverage, minSeqIdentity):
	#out_file = open(idName + '_PDBMatchSummary.out', 'w')
	#print >> out_file,"NPid" + "\t\t\t" + "PDBid" + "\t\t" + "#ofChains" + "\t\t" + "MissingResidues"
	listPDB=[]
	inp = open(NP_id + '.out', 'r')
	for line in inp.readlines():
         if not line.startswith('#'):
		currRow = line.split()
		listPDB.append(currRow)
	pdb_id = ''   # create empty string to hold first PDB id that meets criteria
# 	for i in range(0, len(listPDB)):
#		if (int(listPDB[i][2]) >= minQueryCoverage and int(listPDB[i][3]) >= minSeqIdentity):
#			pdb_id = listPDB[0][0]
#                        print pdb_id
#			actQueryCoverage = listPDB[0][2]
#			actSeqIdentity = listPDB[0][3]
#			matchRun=runPDB(NP_id,pdb_id, actQueryCoverage, actSeqIdentity)
 
	if (int(listPDB[0][2]) >= minQueryCoverage and int(listPDB[0][3]) >= minSeqIdentity):
			pdb_id = listPDB[0][0]
			chain_id = listPDB[0][1]
                        print pdb_id, chain_id
			actQueryCoverage = listPDB[0][2]
			actSeqIdentity = listPDB[0][3]
			seqStartPos = listPDB[0][5]
			seqEndPos = listPDB[0][6]
			sbjStartPos = listPDB[0][7]
			sbjEndPos = listPDB[0][8]
			matchRun=runPDB(NP_id,pdb_id, chain_id,actQueryCoverage, actSeqIdentity, seqStartPos, seqEndPos, sbjStartPos, sbjEndPos)

def matchPDBHighThroughput(NP_id, minQueryCoverage, minSeqIdentity):
	#out_file = open(idName + '_PDBMatchSummary.out', 'w')
	#print >> out_file,"NPid" + "\t\t\t" + "PDBid" + "\t\t" + "#ofChains" + "\t\t" + "MissingResidues"
	listPDB=[]
	inp = open(NP_id + '.out', 'r')
	for line in inp.readlines():
         if not line.startswith('#'):
		currRow = line.split()
		listPDB.append(currRow)
	pdb_id = ''   # create empty string to hold first PDB id that meets criteria
 	for i in range(0, len(listPDB)):
		if (int(listPDB[i][2]) >= minQueryCoverage and int(listPDB[i][3]) >= minSeqIdentity):
			pdb_id = listPDB[i][0]
			actQueryCoverage = listPDB[i][2]
			actSeqIdentity = listPDB[i][3]
			#print(pdb_id)
			break
	if pdb_id.__len__() > 0:
                #http://www.rcsb.org/pdb/rest/getBioAssemblies?structureId=1GAV
		#pdb_output = open(pdb_id + "_sequenceFromPDB.out", "w")
		pdb_file_name = pdb_id + ".pdb"
		pdb_file = open(pdb_file_name, "w")
		url = 'http://www.rcsb.org/pdb/files/%s.pdb' % pdb_id
		pdb = urllib.urlopen(url).read()
		pdb_file.write(pdb)	
		pdb_file.close()

		pdb_file = open(pdb_file_name, "rU")
		structure = PDBParser().get_structure('InPdb', pdb_file)
		model = structure[0]
		numberOfChains = len(model)
		polypeptide = PPBuilder().build_peptides(structure, aa_only = False)
		ppList = []
		if len(polypeptide) == 1:
			ppList.extend(polypeptide[0])
		else:
			for i in range(len(polypeptide)):
				ppList.extend(polypeptide[i])
	
		missingResidues = False
		structureLength = len(ppList)
		for n in range(len(ppList)-1):
			residue_1 = ppList[n]
			residue_id_1 = residue_1.get_id()
			resseq_1 = residue_id_1[1]
			residue_2 = ppList[n+1]
			residue_id_2 = residue_2.get_id()
			resseq_2 = residue_id_2[1]
			if ((resseq_1 + 1) != resseq_2):
				missingResidues = True
				break
				
		print NP_id,pdb_id,missingResidues,numberOfChains
		return (NP_id + "\t" + pdb_id + "\t" + str(numberOfChains) + "\t" + str(missingResidues) + "\t" + str(actQueryCoverage) + "\t" + str(actSeqIdentity) + "\n")



def runPDB(NP_id, pdb_id, chain_id,actQueryCoverage, actSeqIdentity, seqStartPos, seqEndPos, sbjStartPos, sbjEndPos):
	if pdb_id.__len__() > 0:	
		#pdb_output = open(pdb_id + "_sequenceFromPDB.out", "w")
		pdb_file_name = pdb_id + ".pdb"
		if not os.path.exists(pdb_file_name):
			num = 1
                        code = 200
                        url = 'http://www.rcsb.org/pdb/files/'
                        while (code != 404):
                            newUrl = url + pdb_id + ".pdb"+str(num)
                            con = urllib.urlopen(newUrl)
                            code = con.getcode()
                            if(code != 404):
                                print 'There is at least one biological assembly'
                                getPDBFiles(pdb_id)
                                if (num == 1):
                                    os.system('python MakeMultimer.py -c 1 '+ pdb_id + ".pdb")
                                os.rename(pdb_id + '_mm' + str(num) + '.pdb', pdb_id +'_BA_' + str(num) + '.pdb')
                                pdbFileBAName = pdb_id +'_BA_' + str(num) + '.pdb'
				print "I am ready to run readPDBFileforNP"
                                readPDBFileforNP(NP_id, pdb_id, chain_id,actQueryCoverage, actSeqIdentity, seqStartPos, seqEndPos, sbjStartPos, sbjEndPos,pdbFileBAName)
                                
                            elif (code==404) and (num == 1):
                                print 'There are no biological assemblies'    
                                getPDBFiles(pdb_id)
                                pdb_file_name = pdb_id + ".pdb"
                                readPDBFileforNP(NP_id, pdb_id, chain_id,actQueryCoverage, actSeqIdentity, seqStartPos, seqEndPos, sbjStartPos, sbjEndPos,pdb_file_name)
                           
                            num += num + 1
                            
                else:
                        num = 1
                        code = 200
                        url = 'http://www.rcsb.org/pdb/files/'
                        print "url", url
                        while (code != 404):
                            newUrl = url + pdb_id + ".pdb"+str(num)
                            con = urllib.urlopen(newUrl)
                            code = con.getcode()
                            if(code != 404):
                                print 'There is at least one biological assembly'
                                #check if structure has biomt info
                                if (num == 1):
                                    os.system('python MakeMultimer.py -c 1 '+ pdb_id + ".pdb")
                                os.rename(pdb_id + '_mm' + str(num) + '.pdb', pdb_id +'_BA_' + str(num) + '.pdb')
                                pdbFileBAName = pdb_id +'_BA_' + str(num) + '.pdb'
				print "I am ready to run readPDBFileforNP"
                                readPDBFileforNP(NP_id, pdb_id, chain_id,actQueryCoverage, actSeqIdentity, seqStartPos, seqEndPos, sbjStartPos, sbjEndPos,pdbFileBAName)
                                
                            elif (code==404) and (num == 1):
                                print 'There are no biological assemblies'    
                                pdb_file_name = pdb_id + ".pdb"
                                readPDBFileforNP(NP_id, pdb_id, chain_id,actQueryCoverage, actSeqIdentity, seqStartPos, seqEndPos, sbjStartPos, sbjEndPos,pdb_file_name)
                           
                            num += num + 1
#
def runPDBListforBA(NP_id, pdb_id, chain_id):
	if pdb_id.__len__() > 0:	
		#pdb_output = open(pdb_id + "_sequenceFromPDB.out", "w")
		pdb_file_name = pdb_id + ".pdb"
		if not os.path.exists(pdb_file_name):
			num = 1
                        code = 200
                        url = 'http://www.rcsb.org/pdb/files/'
                        while (code != 404):
                            newUrl = url + pdb_id + ".pdb"+str(num)
                            con = urllib.urlopen(newUrl)
                            code = con.getcode()
                            if(code != 404):
                                print 'There is at least one biological assembly'
                                getPDBFiles(pdb_id)
                                if (num == 1):
                                    os.system('python MakeMultimer.py -c 1 '+ pdb_id + ".pdb")
                                os.rename(pdb_id + '_mm' + str(num) + '.pdb', pdb_id +'_BA_' + str(num) + '.pdb')
                                pdbFileBAName = pdb_id +'_BA_' + str(num) + '.pdb'
                                       
                            elif (code==404) and (num == 1):
                                print 'There are no biological assembly'    
                                getPDBFiles(pdb_id)
                                pdb_file_name = pdb_id + ".pdb"

                            num += num + 1   

                else:
                        num = 1
                        code = 200
                        url = 'http://www.rcsb.org/pdb/files/'
                        while (code != 404):
                            newUrl = url + pdb_id + ".pdb"+str(num)
                            con = urllib.urlopen(newUrl)
                            code = con.getcode()
                            if(code != 404):
                                print 'There is at least one biological assembly'
                                #---- check if the structure has BIOMT in it so we need to check because some of the files we have coming from Modeller 
                                if (num == 1):
                                    os.system('python MakeMultimer.py -c 1 '+ pdb_id + ".pdb")
                                os.rename(pdb_id + '_mm' + str(num) + '.pdb', pdb_id +'_BA_' + str(num) + '.pdb')
                                pdbFileBAName = pdb_id +'_BA_' + str(num) + '.pdb'
                                        
                            elif (code==404) and (num == 1):
                                print 'There are no biological assembly'    
                                getPDBFiles(pdb_id)
                                pdb_file_name = pdb_id + ".pdb"

                            num += num + 1

# runPDBList works for pdb file only
def runPDBList(NP_id, pdb_id, chain_id):
	if pdb_id.__len__() > 0:	
		#pdb_output = open(pdb_id + "_sequenceFromPDB.out", "w")
		pdb_file_name = pdb_id + ".pdb"
		if not os.path.exists(pdb_file_name):
			num = 1
                        code = 200
                        url = 'http://www.rcsb.org/pdb/files/'
                        while (code != 404):
                            newUrl = url + pdb_id + ".pdb"+str(num)
                            con = urllib.urlopen(newUrl)
                            code = con.getcode()
                            if(code != 404):
                                print 'There is at least one biological assembly'
                                getPDBFiles(pdb_id)
                                if (num == 1):
                                    os.system('python MakeMultimer.py -c 1 '+ pdb_id + ".pdb")
                                os.rename(pdb_id + '_mm' + str(num) + '.pdb', pdb_id +'_BA_' + str(num) + '.pdb')
                                pdbFileBAName = pdb_id +'_BA_' + str(num) + '.pdb'
                                readPDBFile(NP_id,pdb_id,chain_id,pdbFileBAName)
                                        
                            elif (code==404) and (num == 1):
                                print 'There are no biological assembly'    
                                getPDBFiles(pdb_id)
                                pdb_file_name = pdb_id + ".pdb"
                                readPDBFile(NP_id,pdb_id,chain_id,pdb_file_name)

                            num += num + 1   

                else:
                        num = 1
                        code = 200
                        url = 'http://www.rcsb.org/pdb/files/'
                        while (code != 404):
                            newUrl = url + pdb_id + ".pdb"+str(num)
                            con = urllib.urlopen(newUrl)
                            code = con.getcode()
                            if(code != 404):
                                print 'There is at least one biological assembly'
                                #---- check if the structure has BIOMT in it so we need to check because some of the files we have coming from Modeller 
                                if (num == 1):
                                    os.system('python MakeMultimer.py -c 1 '+ pdb_id + ".pdb")
                                os.rename(pdb_id + '_mm' + str(num) + '.pdb', pdb_id +'_BA_' + str(num) + '.pdb')
                                pdbFileBAName = pdb_id +'_BA_' + str(num) + '.pdb'
                                readPDBFile(NP_id,pdb_id,chain_id,pdbFileBAName)
                                        
                            elif (code==404) and (num == 1):
                                print 'There are no biological assembly'    
                                getPDBFiles(pdb_id)
                                pdb_file_name = pdb_id + ".pdb"
                                readPDBFile(NP_id,pdb_id,chain_id,pdb_file_name)

                            num += num + 1



#------ fetch pdb structure and apply fortran

def readPDBFileforNP(NP_id, pdb_id, chain_id,actQueryCoverage, actSeqIdentity, seqStartPos, seqEndPos, sbjStartPos, sbjEndPos,pdb_file_name):
		pdb_file = open(pdb_file_name, "rU")
		structure = PDBParser().get_structure('InPdb', pdb_file)
		model = structure[0]
		numberOfChains = len(model)
		polypeptide = PPBuilder().build_peptides(structure, aa_only = False)
		ppList = []
		if len(polypeptide) == 1:
			ppList.extend(polypeptide[0])
		else:
			for i in range(len(polypeptide)):
				ppList.extend(polypeptide[i])
	
		missingResidues = False
		structureLength = len(ppList)


#------- sequence-structure id alignment part ---------

## get sequence information #

		SequenceList=[]
		print NP_id
		if not (NP_id == 0):	
			fasta_file_name =NP_id+".fasta"
			inFile=open(fasta_file_name,"rU")
			seqrecs = SeqIO.parse(inFile,'fasta')
			allrecs = list(seqrecs)
			SequenceList.append(allrecs[0].seq)
			#print SequenceList[0]

## get sequence from PDB structure
		        StrSeqList=[]
			for pp in polypeptide:
				StrSeqList.append( pp.get_sequence())
 		               	#print StrSeqList[0]

         		SeqNoListF=[]
          		for i in range(int(seqStartPos),int(seqEndPos)+1):
             			SeqNoListF.append(i)
         		print 'SeqNoList',SeqNoListF


		else:
			print "no fasta file"	
			# I need to put this section because when we give pdb ids. we don't know sequence information. 
			SeqNoList=[]
			i=1
			for i in range (len(ppList)):
				SeqNoList.append(i)
				i+=1	

                print len(ppList)

                SeqNoList =[]
                i=1
                if(SeqNoListF[0] ==1):
                #check structure if it really starts with 1
                       residueInitStr = ppList[0].get_id()[1]
                       residueFinalStr = residueInitStr + len(ppList)
                       print residueInitStr,residueFinalStr
                       if(residueInitStr != SeqNoListF[0]):
                               for i in range(residueInitStr,residueFinalStr+1):
                                       SeqNoList.append(i)
                                       i+=1
                
 
                else:
                        SeqNoList = SeqNoListF
                               
                print 'SeqNoModList',SeqNoList

#------- checking missing residues


		for n in range(len(ppList)-1):
			residue_1 = ppList[n]
			residue_id_1 = residue_1.get_id()
			resseq_1 = residue_id_1[1]
			residue_2 = ppList[n+1]
			residue_id_2 = residue_2.get_id()
			resseq_2 = residue_id_2[1]
			if ((resseq_1 + 1) != resseq_2):
				missingResidues = True
				if(missingResidues == True):#Get pdb files from msv3d
                                   msv3d(NP_id)
				#break

#-----checking structure is complex or not

		if numberOfChains > 1:

#------------------------------------------------
#		# complex
#------------------------------------------------
			Originalcwd = os.getcwd()
			print Originalcwd

#--------------------------------------------------------
#		# calculate dfi for complex structure
#--------------------------------------------------------
			print("complex")
		
			readpdb_file = open(pdb_file_name, "rU")
			writepdb_file = open(pdb_id+'-complex.pdb', "w")
			for line in readpdb_file:
				if line.startswith("ATOM"):
					writepdb_file.write(line)
			writepdb_file.write('END\n')

			writepdb_file.close()
			complexFile=pdb_id+'-complex'
			complexpdbFile= complexFile +'.pdb'
			if not os.path.exists('complexFile'):
				try:os.mkdir(complexFile)
				except:pass
			
			RunFortranFile(complexFile, structureLength)

			SeqNoMonComplexList=[]	
                       #pass first sequence id information
			print numberOfChains
 			for i in range(1):
 			       print 'SeqNoComplexList',i, SeqNoList
			       SeqNoMonComplexList.append(SeqNoList)

                        S2mas = GetDataForComplexStructures(NP_id, complexFile, ppList, SeqNoMonComplexList, numberOfChains)

                        cwd = os.getcwd()
			# get data from DSSP 
			DSSPList=runDSSP(model, complexpdbFile,DSSPPath)

			mergedList = zip(S2mas, DSSPList)



			outname= cwd+"/"+pdb_id+'_dfi-Complex.dat'
			out_file=open(outname,'w')
			out_file.write('Npid \t PDBid \t SeqID \t StrID \t ResName \t dfi \t reldfi \t %dfi \t z-score-dfi \t dsi \t reldsi \t %dsi \t z-score-dsi \t bfactor \t relbfactor \t %bfactor \t ssMotifType \t ASA \t relASA \t %ASA \t phi \t psi \n')
			for i,j in mergedList:
			   out_file.write(str(i[0])+'\t'+str(i[1])+'\t'+str('SeqID')+'\t'+str(i[2])+'\t'+str(i[3])+'\t'+str(i[4])+'\t'+str(i[5])+'\t'+str(i[6])+'\t'+str(i[7])+'\t'+str(i[8])+ \
                               '\t'+ str(i[9])+'\t'+str(i[10])+'\t'+str(i[11])+'\t'+str(i[12])+'\t'+str(i[13])+'\t'+str(i[14])+'\t'+str(j[0])+'\t'+str(j[1])+'\t'+str(j[2])+'\t'+str(j[3])+'\t'+str(j[4])+ \
                               '\t'+str(j[5])+'\n')

    			out_file.close()


#-----------------------------------------------------------------
		# calculate dfi for monomer part of the complex
#-----------------------------------------------------------------
			
			print("monomer from complex")
			try:os.chdir(Originalcwd)
			except:pass
#			ObtainMonomer(NP_id,pdb_id,structure,model,ppList, SeqNoList)
			ObtainMonomerFromComplex(NP_id,pdb_id, chain_id,structure,model,ppList,SeqNoMonComplexList)
			
			os.chdir(Originalcwd)

		else:
#----------------------------------------------------------------
#		# calculate dfi for monomer
#----------------------------------------------------------------

			cwd = os.getcwd()
			print cwd

			print("monomer")
#			try:os.mkdir(pdb_id)
#			except:pass

#			print numberOfChains
# 			for i in range(1):
# 			       print 'SeqNoComplexList',i, SeqNoList
#			       SeqNoMonComplexList.append(SeqNoList)

        		no_of_models = len(structure.get_list())
			print no_of_models
			if(no_of_models>1):
				print 'this is NMR'
				ObtainMonomerNMR(NP_id,pdb_id,structure,model,ppList, SeqNoList)
			else:
				ObtainMonomer(NP_id,pdb_id,structure,model,ppList, SeqNoList)
			
			try:os.chdir(cwd)
			except:pass



		#writemas(S2mas)
		return (str(NP_id) + "\t" + str(pdb_id) + "\t" + str(numberOfChains) + "\t" + str(missingResidues) + "\t" + str(actQueryCoverage) + "\t" + str(actSeqIdentity) + "\n") 
		
		# out_file.write(idName + "\t\t")
		# out_file.write(pdb_id + "\t\t")
		# out_file.write(str(numberOfChains) + "\t\t\t")
		# out_file.write(str(missingResidues) + "\n")
					
	# else:
		# return (NP_id + "\t\tCoverage and identity minima not met\n")
		#out_file.write('Coverage and identity minima not met')
	#out_file.close()
	
#matchPDB("NP_000006", 90, 90)




def readPDBFile(NP_id,pdb_id,chain_id,pdb_file_name):
		pdb_file = open(pdb_file_name, "rU")
		structure = PDBParser().get_structure('InPdb', pdb_file)
		model = structure[0]
		numberOfChains = len(model)
		polypeptide = PPBuilder().build_peptides(structure, aa_only = False)
		ppList = []
		if len(polypeptide) == 1:
			ppList.extend(polypeptide[0])
		else:
			for i in range(len(polypeptide)):
				ppList.extend(polypeptide[i])
	
		missingResidues = False
		structureLength = len(ppList)
                ##fixing##
		SequenceList=[]
		print NP_id
		if (NP_id == 0):	
			print "no fasta file"	
			# I need to put this section because when we give pdb ids. we don't know sequence information. 
			SeqNoList=[]
			j=1
			#designed to read only the residues in atom,  will leave out water etc
                        for model in structure:
                         for chain in model:
                          for residue in chain:
                                   temp = residue.get_id()
                                   if temp[0] == ' ':
                                           SeqNoList.append(j)
                                           j+=1
                        print (j-1)
                      
                                   
                        
			#print ppList
			#for i in range (len(ppList)):
			#	SeqNoList.append(i)
			#	i+=1
				 


#------- checking missing residues

		for n in range(len(ppList)-1):
			residue_1 = ppList[n]
			residue_id_1 = residue_1.get_id()
			resseq_1 = residue_id_1[1]
			residue_2 = ppList[n+1]
			residue_id_2 = residue_2.get_id()
			resseq_2 = residue_id_2[1]
			if ((resseq_1 + 1) != resseq_2):
				missingResidues = True
				break
               
#-----checking structure is complex or not

		if numberOfChains > 1:

#------------------------------------------------
#		# complex
#------------------------------------------------
			Originalcwd = os.getcwd()
			print Originalcwd

#--------------------------------------------------------
#		# calculate dfi for complex structure
#--------------------------------------------------------
			print("complex")
		
			readpdb_file = open(pdb_file_name, "rU")
			writepdb_file = open(pdb_id+'-complex.pdb', "w")
			for line in readpdb_file:
				if line.startswith("ATOM"):
					writepdb_file.write(line)
			writepdb_file.write('END\n')

			writepdb_file.close()
			complexFile=pdb_id+'-complex'
			complexpdbFile= complexFile +'.pdb'
			if not os.path.exists('complexFile'):
				try:os.mkdir(complexFile)
				except:pass
			
			RunFortranFile(complexFile, structureLength)

			SeqNoComplexList=[]
			
			for i in range(numberOfChains):
				SeqNoComplexList.append(SeqNoList)

	        	S2mas = GetDataForComplexStructures(NP_id, complexFile, ppList, SeqNoComplexList, numberOfChains)

                        cwd = os.getcwd()
			# get data from DSSP 
			DSSPList=runDSSP(model, complexpdbFile,DSSPPath)

			mergedList = zip(S2mas, DSSPList)	

			outname= cwd+"/"+pdb_id+'_dfi-Complex.dat'
			out_file=open(outname,'w')
			out_file.write('Npid \t PDBid \t SeqID \t StrID \t ResName \t dfi \t reldfi \t %dfi \t z-score-dfi \t dsi \t reldsi \t %dsi \t z-score-dsi \t bfactor \t relbfactor \t %bfactor \t ssMotifType \t ASA \t relASA \t %ASA \t phi \t psi \n')
			for i,j in mergedList:
			   out_file.write(str(i[0])+'\t'+str(i[1])+'\t'+str('SeqID')+'\t'+str(i[2])+'\t'+str(i[3])+'\t'+str(i[4])+'\t'+str(i[5])+'\t'+str(i[6])+'\t'+str(i[7])+'\t'+str(i[8])+ \
                               '\t'+ str(i[9])+'\t'+str(i[10])+'\t'+str(i[11])+'\t'+str(i[12])+'\t'+str(i[13])+'\t'+str(i[14])+'\t'+str(j[0])+'\t'+str(j[1])+'\t'+str(j[2])+'\t'+str(j[3])+'\t'+str(j[4])+ \
                               '\t'+str(j[5])+'\n')

      			out_file.close()

		
#---------------------------------------------------------------
		# calculate dfi for monomer part of the complex
#--------------------------------------------------------------
			
			print("monomer from complex")
			try:os.chdir(Originalcwd)
			except:pass

#			ObtainMonomer(NP_id,pdb_id,structure,model,ppList, SeqNoList)

  			ObtainMonomerFromComplex(NP_id,pdb_id, chain_id,structure,model,ppList,SeqNoComplexList)
			
			os.chdir(Originalcwd)

		else:
#------------------------------------------------
#		# calculate dfi for monomer
#------------------------------------------------

			cwd = os.getcwd()
			print cwd

			print("monomer")
#			try:os.mkdir(pdb_id)
#			except:pass

        		no_of_models = len(structure.get_list())
			print no_of_models
			if(no_of_models>1):
				print 'this is NMR'
				ObtainMonomerNMR(NP_id,pdb_id,structure,model,ppList, SeqNoList)
			else:
				ObtainMonomer(NP_id,pdb_id,structure,model,ppList, SeqNoList)
			
			try:os.chdir(cwd)
			except:pass



		#writemas(S2mas)
		return (str(NP_id) + "\t" + str(pdb_id) + "\t" + str(numberOfChains) + "\t" + str(missingResidues) +  "\n") 
		
		# out_file.write(idName + "\t\t")
		# out_file.write(pdb_id + "\t\t")
		# out_file.write(str(numberOfChains) + "\t\t\t")
		# out_file.write(str(missingResidues) + "\n")
					
	# else:
		# return (NP_id + "\t\tCoverage and identity minima not met\n")
		#out_file.write('Coverage and identity minima not met')
	#out_file.close()
	
#matchPDB("NP_000006", 90, 90)








#----------- obtain Monomer structure from NMR

def ObtainMonomerNMR(NP_id,pdb_id,structure,model,ppList, SeqNoList):
		chainID_list = []
		for mod in model:
			chainID_list.append(mod.get_id())
		chainID = chainID_list[0]

 		pdbFile = pdb_id + '.pdb'
       		f = open(pdbFile,'rb')
        	pdb = f.readlines()
        	f.close()
                
                cwd = os.getcwd()
                shutil.rmtree(cwd +'/'+pdb_id+'-monomer', True)
                
		models = splitNMR(pdb)
		monomer_pdb = pdb_id + '_monomer'
		monomer_pdbFile =monomer_pdb+'.pdb'
		g= open(monomer_pdbFile,"w")
		g.writelines(models[0])
		g.close
		g.flush()



		readpdb_file = open(monomer_pdbFile, "rb")
		writepdb_file = open(pdb_id+'-monomer.pdb', "w")
		for line in readpdb_file.readlines():
			if line.startswith("ATOM"):
				writepdb_file.write(line)
		writepdb_file.write('END\n')
		writepdb_file.close()

		monomerFile=pdb_id+'-monomer'
		monomerpdbFile= monomerFile+'.pdb'
		if not os.path.exists('monomerFile'):
			try:os.mkdir(monomerFile)
			except:pass
                os.remove(monomer_pdbFile)
		
                structureLength = len(ppList)


##-- to do: get binding and catalytic data 
#		monomerpdbsumList = pdbsum(pdb_id,chainID, ppList)


		RunFortranFile(monomerFile, structureLength)
	        S2mas = GetDataForNMRStructures(NP_id, monomerFile, ppList, SeqNoList)

                cwd = os.getcwd()

		# get data from DSSP 
		DSSPList=runDSSP(model, monomerpdbFile,DSSPPath)

		mergedList = zip(S2mas, DSSPList)	

#		os.chdir(cwd)
		outname= cwd+"/"+ monomer_pdb +'_dfi.dat'
		out_file=open(outname,'w')
		out_file.write('Npid \t PDBid \t SeqID \t StrID \t ResName \t dfi \t reldfi \t %dfi \t z-score dfi \t dsi \t reldsi \t %dsi \t z-score dsi \t bfactor \t relbfactor \t %bfactor \t ssMotifType \t ASA \t relASA \t %ASA \t phi \t psi \n')    

		for i,j in mergedList:
                         out_file.write(str(i[0])+'\t'+str(i[1])+'\t'+str(i[2])+'\t'+str(i[3])+'\t'+str(i[4])+'\t'+str(i[5])+'\t'+str(i[6])+'\t'+str(i[7])+ '\t' + str(i[8]) + \
                                '\t'+ str(i[9]) + '\t' + str(i[10]) + '\t' + str(i[11]) + '\t'+ str(i[12]) + '\t' + str('0') + '\t'+ str('0') +'\t'+ str('0') +'\t'+str(j[0])+'\t'+str(j[1])+ \
                                '\t'+str(j[2])+'\t'+str(j[3])+'\t'+str(j[4])+'\t'+str(j[5])+'\t'+'\n')
         	out_file.close()


                

def splitNMR(pdb):
    """
    Split each model in an NMR pdb file into its own pdb.
    """

    to_strip = ["ENDMDL","MASTER"]
    pdb = [l for l in pdb if l[0:6] not in to_strip]
    pdb_hash = [(l[0:6],i) for i, l in enumerate(pdb)]
    pdb_hash = [x[1] for x in pdb_hash if x[0] == "MODEL "]
    pdb_hash.append(len(pdb))

    all_models = []
    for i in range(1,len(pdb_hash)):
        all_models.append(pdb[pdb_hash[i-1]:pdb_hash[i]])
       
    return all_models 



#----------- obtain Monomer structure from X-ray structure

def ObtainMonomer(NP_id,pdb_id,structure,model,ppList, SeqNoList):
                
		chainID_list = []
		for mod in model:
			chainID_list.append(mod.get_id())
		chainID = chainID_list[0]
# this part is kind of awkward 

		monomer_pdb = pdb_id + '_monomer'
		monomer_pdbFile =monomer_pdb+'.pdb'
		io=PDBIO()
		io.set_structure(structure)
		io.save(monomer_pdbFile)			
#		io.save(monomer_pdbFile, select=NotDisordered(), chain_select())			
		

		readpdb_file = open(monomer_pdbFile, "rU")
		writepdb_file = open(pdb_id+'-monomer.pdb', "w")
		for line in readpdb_file:
			if line.startswith("ATOM"):
				writepdb_file.write(line)
		writepdb_file.write('END\n')

		writepdb_file.close()
		readpdb_file.close()
		monomerFile=pdb_id+'-monomer'
		monomerpdbFile= monomerFile+'.pdb'
		if not os.path.exists('monomerFile'):
			try:os.mkdir(monomerFile)
			except:pass
                os.remove(monomer_pdbFile)			
		
#		structure = PDBParser().get_structure('InPdb', monomerpdbFile)
#		model = structure[0]
#		polypeptide = PPBuilder().build_peptides(structure, aa_only = False)
#		MonomerppList = []
#		if len(polypeptide) == 1:
#			MonomerppList.extend(polypeptide[0])
		

                structureLength = len(ppList)
                sequenceLength = len(SeqNoList)

		print 'length of polypeptide in structure', structureLength
 		print 'length of polypeptide in sequence', sequenceLength

# 		if(structureLength != sequenceLength):
#                        structureLength = sequenceLength
               
##-- to do: get binding and catalytic data 
		#monomerpdbsumList = pdbsum(pdb_id,chainID, ppList)
		#print monomerpdbsumList
                #bindResList=[]

                #j=0
                #for i in range(structureLength):
                #        print ppList[i].get_id()[1]
                #        if (ppList[i].get_id()[1]== monomerpdbsumList[j]):
                #                bindResList.append(1)
                #                j=j+1
                #        else:
                #                bindResList.append(0)
                #                j=j=1
                #print bindResList        
                        


		RunFortranFile(monomerFile, structureLength)
	        S2mas = GetData(NP_id, monomerFile, ppList, SeqNoList)

		# get data from DSSP
                cwd = os.getcwd()
		DSSPList=runDSSP(model, monomerpdbFile,DSSPPath)

		print "DSSPList", DSSPList

#		mergedList = zip(S2mas, DSSPList, monomerpdbsumList)	
		mergedList = zip(S2mas, DSSPList)
                
#                currentPath=os.getcwd()
#		print currentPath
#		os.chdir(cwd)
		outname= cwd+"/"+ monomer_pdb +'_dfi.dat'
		out_file=open(outname,'w')


		out_file.write('Npid \t PDBid \t SeqID \t StrID \t ResName \t dfi \t reldfi \t %dfi \t z-score-dfi \t dsi \t reldsi \t %dsi \t z-score-dsi \t bfactor \t relbfactor \t %bfactor \t ssMotifType \t ASA \t relASA \t %ASA \t phi \t psi \n')    
		for i,j in mergedList:
                         out_file.write(str(i[0])+'\t'+str(i[1])+'\t'+str(i[2])+'\t'+str(i[3])+'\t'+str(i[4])+'\t'+str(i[5])+'\t'+str(i[6])+'\t'+str(i[7])+ '\t' + str(i[8]) + \
                                '\t'+ str(i[9]) + '\t' + str(i[10]) + '\t' + str(i[11]) + '\t'+ str(i[12]) + '\t' + str(i[13]) + '\t'+ str(i[14]) +'\t'+ str(i[15]) +'\t'+str(j[0])+'\t'+str(j[1])+ \
                                '\t'+str(j[2])+'\t'+str(j[3])+'\t'+str(j[4])+'\t'+str(j[5])+'\n')
                         
         	out_file.close()

                
                

 


def ObtainMonomerFromComplex(NP_id,pdb_id,chain_id,structure,model,ppList,SeqNoList):

		print "i am chain:", chain_id
# this part is kind of awkward 

		monomer_pdb = pdb_id + '_monomer'
		monomer_pdbFile =monomer_pdb+'.pdb'
		io=PDBIO()
		io.set_structure(structure)
#		io.save(monomer_pdbFile, chain_select())			
		io.save(monomer_pdbFile, select=chain_select())			
		

		readpdb_file = open(monomer_pdbFile, "rU")
		writepdb_file = open(pdb_id+'-monomer.pdb', "w")
		for line in readpdb_file:
			if line.startswith("ATOM"):
				writepdb_file.write(line)
		writepdb_file.write('END\n')

		writepdb_file.close()
		monomerFile=pdb_id+'-monomer'
		monomerpdbFile= monomerFile+'.pdb'
		if not os.path.exists('monomerFile'):
			try:os.mkdir(monomerFile)
			except:pass
		
		os.remove(monomer_pdbFile)

 
		structure = PDBParser().get_structure('InPdb', monomerpdbFile)
		polypeptide = PPBuilder().build_peptides(structure[0], aa_only = False)
		print polypeptide[0]
		for model in structure:
		  for chain in model:
		   i = 0
		   for residue in chain:
    			i += 1
		print i

#		model = structure[0]
#		polypeptide = PPBuilder().build_peptides(structure, aa_only = False)
		MonomerppList = []
		MonomerppList.extend(polypeptide[0])
#
#		if len(polypeptide) == 1:
#			MonomerppList.extend(polypeptide[0])
#		else:
#			print polypeptide[0]
#			MonomerppList.extend(polypeptide[0])
		
		chainID = chain_id

	        structureLength = len(MonomerppList)
		print 'length of polypeptide in monomer', structureLength

 		SeqNoList=[]
		i=1
		for i in range (len(MonomerppList)):
				SeqNoList.append(i)
				i+=1
		print 'checking again',SeqNoList

#		structureLength = i

		print structureLength

#		monomerpdbsumList = pdbsum(pdb_id,chainID, MonomerppList)

#		print monomerpdbsumList

		RunFortranFile(monomerFile, structureLength)
	        S2mas = GetData(NP_id, monomerFile, MonomerppList, SeqNoList)

		cwd = os.getcwd()
		# get data from DSSP 
		DSSPList=runDSSP(model, monomerpdbFile,DSSPPath)

		mergedList = zip(S2mas, DSSPList)	
#		mergedList = zip(S2mas, DSSPList, SeqNoList)	


#		os.chdir(cwd)
		outname= cwd+"/"+ monomer_pdb +'_dfi.dat'
		out_file=open(outname,'w')
		out_file.write('Npid \t PDBid \t SeqID \t StrID \t ResName \t dfi \t reldfi \t %dfi \t z-score dfi \t dsi \t reldsi \t %dsi \t z-score dsi \t bfactor \t relbfactor \t %bfactor \t ssMotifType \t ASA \t relASA \t %ASA \t phi \t psi \n')    

		for i,j in mergedList:
	 	 out_file.write(str(i[0])+'\t'+str(i[1])+'\t'+str(i[2])+'\t'+str(i[3])+'\t'+str(i[4])+'\t'+str(i[5])+'\t'+str(i[6])+'\t'+str(i[7])+ '\t' + str(i[8]) + \
                                '\t'+ str(i[9]) + '\t' + str(i[10]) + '\t' + str(i[11]) + '\t'+ str(i[12]) + '\t' + str(i[13]) + '\t'+ str(i[14]) +'\t'+ str(i[15]) +'\t'+str(j[0])+'\t'+str(j[1])+ \
                                '\t'+str(j[2])+'\t'+str(j[3])+'\t'+str(j[4])+'\t'+str(j[5])+'\n')
         	out_file.close()



#------ creates specific fortran file with the input file 

def RunFortranFile(pdb_id, numResidues):
	File = pdb_id + '.pdb' 
	print File, pdb_id,numResidues

	cwd = os.getcwd()
	checkPath = cwd+'/'+pdb_id
	print checkPath
	#make folder to move everything into, if folder exists then just move things into that folder

                
# write a file for CA atoms
        numResiduesCA = createCAfile(pdb_id)
        fileCA = pdb_id + '-CA.pdb'
        shutil.move(fileCA,pdb_id)
        
	if not os.path.exists(checkPath+'/'+File):
		shutil.move(File, pdb_id)

	CAFile_name = pdb_id +'-CA'

	

        # check if there is already created fortran files
	input_name = pdb_id + 'input.parms'

	# first check these files are already in the pdb_id folder
#        os.chdir(pdb_id)
#	if not os.path.exists(checkPath+'/'+input_name):
	input = open(input_name, 'w')
	fort = pdb_id + 'prs.f'
	destination = open (fort, 'w')
	print str(numResiduesCA)
	input.write('1\t' + str(numResiduesCA) + '\t1\t' + str(numResiduesCA) + '\t1')
	input.close()
 
	#open fortran file and make a copy, then rewrite it and copy over the original
	source = open('prs.f', 'r')
	
	i=1
	linenums=[25,170,182,183,191,193,194,195,204,205,206,207,209,210]
	for line in source:
		if i not in linenums:
			destination.write(line[:-2]+'\n')
		else:
			if i==25:
				destination.write( line[:15]+str(numResiduesCA)+line[18:-2] + "\n" )
				linenums.remove(i)
			elif i==170:
				destination.write( line[:15] + input_name + "')" + '\n')
				linenums.remove(i)
			elif i==182:
				destination.write( line[:15] + CAFile_name + line[17:-2] + "\n")
				linenums.remove(i)
			elif i==183:
				destination.write( line[:15] + CAFile_name + line[17:-2] + "\n")
				linenums.remove(i)
			else:
				destination.write( line[:15] + pdb_id + line[17:-2] + "\n")
				linenums.remove(i)
		i+=1
	source.close()
	destination.close()
	os.system('chmod u+x '+fort)
#	print os.path.exists(pdb_id)
#	if not os.path.exists(pdb_id):
#		if os.path.isdir(pdb_id):
#	if os.path.exists(pdb_id):
#		os.chdir(pdb_id)
#		#in the directory with the pdb file, the input parms, and the fortran file, now compile and run
#		os.system('g77 -o run '+fort)
#		os.system('./run')
#		#move back to parent directory
#		os.chdir('..')
#	else:
	shutil.move(input_name, pdb_id)
	shutil.move(fort,pdb_id)
	os.chdir(pdb_id)
	#in the directory with the pdb file, the input parms, and the fortran file, now compile and run
	os.system('g77 -o run '+fort)
	os.system('./run')
	#move back to parent directory
	os.chdir('..')


#------- Data analyze part

def GetData(npid, pdbid, pp_prot, SeqNoList):
	File=pdbid+'.pdb'
	print File
	#call getrightPDB.py have it return PDB with correct numbering of residues
	#os.system('getrightPDB.py '+pdbid)
	#get res list to cycle through and add to dict later
	# prot=PPeptide()
	# prot.readPDB(pdbid)
	# prot.renumber()
        os.chdir(pdbid)

        s2z_score = []
	s2l,s2rel,s2perc=[],[],[]    #make dict with each residue num and corresponding S2 value
	bfactorl,bfactorrel,bfactorperc=[],[],[]    #make dict with each residue num and corresponding factor value

	s2name=pdbid+'-S2-Avg.dat'
	s2file=open(s2name,'r')
	for line in s2file:
		s2l.append(float(line))
	s2file.close()

        # checks for missing residue
#        if len(pp_prot) < len(SeqNoList):
#                print("ERROR:  There is a missing residue")

	#calc rel S2's and %rank S2's and add to dict
	#get avg
	
        S2avg=sum(s2l)/len(s2l)

	###gets percentage for each element of s2###
	for m in s2l:
                count=0
		s2rel.append(m/S2avg)
		for y in s2l:
			if y<=m:
				count+=1
		s2perc.append(float(count)/float(len(s2l)))

        ##calculated standard deviation
	sum_diff = 0
	for m in s2l:
                diff = (m - S2avg)
                diff_sqr = diff * diff
                sum_diff += diff_sqr

        stndrdDev = sqrt(sum_diff/len(s2l))        

        ##calculates the z score
	for m in s2l:
                s2z_score.append(float(m-S2avg)/float(stndrdDev))
  
        #####data for S1#####
        s1z_score = []
	s1l,s1rel,s1perc=[],[],[]    #make dict with each residue num and corresponding S2 value

	s1name=pdbid+'-S1-Avg.dat'
	s1file=open(s1name,'r')
	for line in s1file:
		s1l.append(float(line))
	s1file.close()
	
	#calc rel S1's and %rank S1's and add to dict
	###get avg for s1###
	S1avg=sum(s1l)/len(s1l)

	###gets percentage for each element of s1###
	for m in s1l:
		count=0
		s1rel.append(m/S1avg)
		for y in s1l:
			if y<=m:
				count+=1
		s1perc.append(float(count)/float(len(s1l)))

        ##calculated standard deviation
	sum_diff = 0
	for m in s1l:
                diff = (m - S1avg)
                diff_sqr = diff * diff
                sum_diff += diff_sqr

        stndrdDev = sqrt(sum_diff/len(s1l))        

        ##calculates the z score
	for m in s1l:
                s1z_score.append(float(m-S1avg)/float(stndrdDev))
              

	p=PDBParser()

# get bfactor values for each residue and calculate percentile
	structure=p.get_structure('InPdb', File) 
	for model in structure:
	 for chain in model:
	  for residue in chain:
	   for atom in residue:
	     if(atom==residue['CA']):
		bfactorl.append(float(atom.get_bfactor()))
	#get avg
	print "check bfactor",bfactorl
        bfactoravg =sum(bfactorl)/len(bfactorl)
	print "check bfactorAVG",bfactoravg
	if(bfactoravg==0.0):
           for m in bfactorl:
		bfactorrel.append(0.0)
                bfactorperc.append(0.0)
           print "check bfactorrel", bfactorrel
           print "check bfactorperc", bfactorperc
        else:	   
           for m in bfactorl:
		count=0
		bfactorrel.append(m/bfactoravg)
		for y in bfactorl:
			if y<=m:
				count+=1
		bfactorperc.append(float(count)/float(len(bfactorl)))

        print len(bfactorrel),len(bfactorperc)

	seqID=[]
	i=0
 	for mm in s2l:
	     seqID.append(SeqNoList[i])
	     i+=1

#        print "length", float(len(bfactoravg)),float(len(bfactorrel)),float(len(bfactorperc))


 	
	i=0
	S2mas=[]
	for n in s2l:
#		S2mas.append([npid,pdbid, seqID[i],pp_prot[i].get_id()[1],pp_prot[i].get_resname(),n,s2rel[i],s2perc[i]])
#		S2mas.append([pp_prot.res_list[i][1],n,s2rel[i],s2perc[i]])
#		S2mas.append([npid,pdbid,seqID[i],pp_prot[i].get_id()[1],pp_prot[i].get_resname(),n,s2rel[i],s2perc[i],s2z_score[i],s1l[i], s1rel[i],s1perc[i],s1z_score[i],bfactorl[i],bfactorrel[i],bfactorperc[i]])
		S2mas.append([npid,pdbid])

		i+=1
	#before returning data let's write it out- There is a problem here

#        replaceBValues(File, pdbid, s2perc)

	#now return list of values 
	return S2mas

#- Data analyze for specific NMR structure

def GetDataForNMRStructures(npid, pdbid, pp_prot, SeqNoList):
	File=pdbid+'.pdb'
	print File
        os.chdir(pdbid)

        #####data for S2#####
	s2z_score = []
	s2l,s2rel,s2perc=[],[],[]    #make dict with each residue num and corresponding S2 value
#	bfactorl,bfactorrel,bfactorperc=[],[],[]    #make dict with each residue num and corresponding factor value

	s2name=pdbid+'-S2-Avg.dat'
	s2file=open(s2name,'r')
	for line in s2file:
		s2l.append(float(line))
	s2file.close()
	
	#calc rel S2's and %rank S2's and add to dict
	###get avg for s2###
	S2avg=sum(s2l)/len(s2l)

	###gets percentage for each element of s2###
	for m in s2l:
		count=0
		s2rel.append(m/S2avg)
		for y in s2l:
			if y<=m:
				count+=1
		s2perc.append(float(count)/float(len(s2l)))

        ##calculated standard deviation
	sum_diff = 0
	for m in s2l:
                diff = (m - S2avg)
                diff_sqr = diff * diff
                sum_diff += diff_sqr

        stndrdDev = sqrt(sum_diff/len(s2l))        

        ##calculates the z score
	for m in s2l:
                s2z_score.append(float(m-S2avg)/float(stndrdDev))

        #####data for S1#####
        s1z_score = []
	s1l,s1rel,s1perc=[],[],[]    #make dict with each residue num and corresponding S2 value

	s1name=pdbid+'-S1-Avg.dat'
	s1file=open(s1name,'r')
	for line in s1file:
		s1l.append(float(line))
	s1file.close()
	
	#calc rel S1's and %rank S1's and add to dict
	###get avg for s1###
	S1avg=sum(s1l)/len(s1l)

	###gets percentage for each element of s1###
	for m in s1l:
		count=0
		s1rel.append(m/S1avg)
		for y in s1l:
			if y<=m:
				count+=1
		s1perc.append(float(count)/float(len(s1l)))

        ##calculated standard deviation
	sum_diff = 0
	for m in s1l:
                diff = (m - S1avg)
                diff_sqr = diff * diff
                sum_diff += diff_sqr

        stndrdDev = sqrt(sum_diff/len(s1l))        

        ##calculates the z score
	for m in s1l:
                s1z_score.append(float(m-S1avg)/float(stndrdDev))

                
	p=PDBParser()

	

# get bfactor values for each residue and calculate percentile
#	structure=p.get_structure('InPdb', File) 
#	for model in structure:
#	 for chain in model:
#	  for residue in chain:
#	   for atom in residue:
#	     if(atom==residue['CA']):
#		bfactorl.append(float(atom.get_bfactor()))
	#get avg
#	bfactoravg =sum(bfactorl)/len(bfactorl)
#	for m in bfactorl:
#		count=0
#		bfactorrel.append(m/bfactoravg)
#		for y in bfactorl:
#			if y<=m:
#				count+=1
#		bfactorperc.append(float(count)/float(len(bfactorl)))
#        replaceBValues(File, pdbid, s2perc)

#	seqID=[]
#	i=0
#	for mm in s2l:
#	     seqID.append(SeqNoList[i])
#	     i+=1
        print 'SeqNoLisinNMR',SeqNoList

	i=0
	S2mas=[]
	for n in s2l:
		S2mas.append([npid,pdbid, 'seqID[i]',pp_prot[i].get_id()[1],pp_prot[i].get_resname(),n,s2rel[i],s2perc[i],s2z_score[i],s1l[i], s1rel[i],s1perc[i],s1z_score[i]])
#		S2mas.append([pp_prot.res_list[i][1],n,s2rel[i],s2perc[i]])
#		S2mas.append([npid,pdbid,SeqNoList,pp_prot[i].get_id()[1],pp_prot[i].get_resname(),n,s2rel[i],s2perc[i],bfactorl[i],bfactorrel[i],bfactorperc[i]])
		i+=1
	#before returning data let's write it out

	#now return list of values 
	return S2mas


def GetDataForComplexStructures(npid, pdbid, pp_prot, SeqNoComplexList, numberOfChains):
	File=pdbid+'.pdb'
	print File
        os.chdir(pdbid)
	
	s2z_score = []
	s2l,s2rel,s2perc=[],[],[]    #make dict with each residue num and corresponding S2 value
	bfactorl,bfactorrel,bfactorperc=[],[],[]  

	s2name=pdbid+'-S2-Avg.dat'
	s2file=open(s2name,'r')
	for line in s2file:
		s2l.append(float(line))
	s2file.close()

	#calc rel S2's and %rank S2's and add to dict
	#get avg
	S2avg=sum(s2l)/len(s2l)

	#gets percentage for each element of s2
	for m in s2l:
		count=0
		s2rel.append(m/S2avg)
		for y in s2l:
			if y<=m:
				count+=1
		s2perc.append(float(count)/float(len(s2l)))

        ##calculated standard deviation
	sum_diff = 0
	for m in s2l:
                diff = (m - S2avg)
                diff_sqr = diff * diff
                sum_diff += diff_sqr

        stndrdDev = sqrt(sum_diff/len(s2l))        

        ##calculates the z score
	for m in s2l:
                s2z_score.append(float(m-S2avg)/float(stndrdDev))

        #####data for S1#####
        s1z_score = []
	s1l,s1rel,s1perc=[],[],[]    #make dict with each residue num and corresponding S2 value

	s1name=pdbid+'-S1-Avg.dat'
	s1file=open(s1name,'r')
	for line in s1file:
		s1l.append(float(line))
	s1file.close()
	
	#calc rel S1's and %rank S1's and add to dict
	###get avg for s1###
	S1avg=sum(s1l)/len(s1l)

	###gets percentage for each element of s1###
	for m in s1l:
		count=0
		s1rel.append(m/S1avg)
		for y in s1l:
			if y<=m:
				count+=1
		s1perc.append(float(count)/float(len(s1l)))

        ##calculated standard deviation
	sum_diff = 0
	for m in s1l:
                diff = (m - S1avg)
                diff_sqr = diff * diff
                sum_diff += diff_sqr

        stndrdDev = sqrt(sum_diff/len(s1l))        

        ##calculates the z score
	for m in s1l:
                s1z_score.append(float(m-S1avg)/float(stndrdDev))



	p=PDBParser()


#	seqID=[]
#	for k in range(numberOfChains):
#	 i=0
#	 for mm in s2l:
#	     print i, SeqNoComplexList[k][i]
#	     seqID.append(SeqNoComplexList[k][i])
#	     i+=1	

# get bfactor values for each residue and calculate percentile
	structure=p.get_structure('InPdb', File) 
	for model in structure:
	 for chain in model:
	  for residue in chain:
	   for atom in residue:
	     if(atom==residue['CA']):
		bfactorl.append(float(atom.get_bfactor()))
	#get avg
	bfactoravg =sum(bfactorl)/len(bfactorl)
	for m in bfactorl:
		count=0
		bfactorrel.append(m/bfactoravg)
		for y in bfactorl:
			if y<=m:
				count+=1
		bfactorperc.append(float(count)/float(len(bfactorl)))

#        replaceBValues(File, pdbid, s2perc)


	i=0
	S2mas=[]
	for n in s2l:
		S2mas.append([npid,pdbid,pp_prot[i].get_id()[1],pp_prot[i].get_resname(),n,s2rel[i],s2perc[i],s2z_score[i],s1l[i], s1rel[i],s1perc[i],s1z_score[i],bfactorl[i],bfactorrel[i],bfactorperc[i]])


#		S2mas.append([pp_prot.res_list[i][1],n,s2rel[i],s2perc[i]])
#		S2mas.append([npid,pdbid,SeqNoList,pp_prot[i].get_id()[1],pp_prot[i].get_resname(),n,s2rel[i],s2perc[i],bfactorl[i],bfactorrel[i],bfactorperc[i]])
		i+=1
	#before returning data let's write it out

	#now return list of values 
	return S2mas


	
def writemas(S2Dict):
	#write out all data to a single file
	MasterOut=open('S2-Compiled-Data.dat','w')
	#sort dict according to num of residues in each structure
	SortedS2List=sorted(S2Dict.items(), key=lambda x: len(x[1][0]), reverse=True)
	#output sorted dict into a single master file
	i=0
	while i<len(SortedS2List):
		MasterOut.write("%4s RESI    S2  REL_S2  PRANK_S2  "%(SortedS2List[i][0]))
		i+=1
	MasterOut.write('\n')
	j,k=0,0
	while k<len(SortedS2List[0][1]):
		j=0
		while j<len(SortedS2List):
			try: MasterOut.write("     %4d %5.4f %5.4f  %5.4f    "%(SortedS2List[j][1][k][0],SortedS2List[j][1][k][1],SortedS2List[j][1][k][2],SortedS2List[j][1][k][3]))
			except: break
			j+=1
		MasterOut.write('\n')
		k+=1
	MasterOut.close()


#--- run DSSP module to get ASA and secondary structural information

def runDSSP(model,File,DSSPPath):
		if not DSSPPath:
		  print 'DSSP path need to be set'
		cwd=os.getcwd()
		File1 = cwd+'/'+File
		d = DSSP(model, File1, DSSPPath)
		ssMotifType,solventAcc,relSolventAcc,solventAccperc,phi,psi =[],[],[],[],[],[]

		for key in d.keys(): 
          	 ssMotifType.append(str(d[key][1]))
		 solventAcc.append(float(d[key][2]))
          	 relSolventAcc.append(float(d[key][3]))
          	 phi.append(float(d[key][4]))
          	 psi.append(float(d[key][5]))
	  
		for m in solventAcc:
			count=0
			for y in solventAcc:
				if y<=m:
					count+=1
			solventAccperc.append(float(count)/float(len(solventAcc)))
#			print solventAcc,solventAccperc


		ssMotifName=[]
		# annotate DSSPmotif type in a right way
		for letter in ssMotifType:
	 		if(letter=='H' or letter=='G' or letter=='I'):
				ssMotifName.append('Helix')
	 		if(letter=='E' or letter=='B'):
	        		ssMotifName.append('Strand')	
	 		if(letter=='T' or letter=='S' or letter=='-'):
				ssMotifName.append('Loop')	


 
		i=0
		DSSPList=[]
		for n in solventAcc:
			DSSPList.append([ ssMotifName[i],solventAcc[i], relSolventAcc[i],solventAccperc[i],phi[i],psi[i]])
			i+=1
		return DSSPList


def pdbsum(pdb_id,chainName,pp_prot):
#------------------------------------------------
## Get binding data from PDBSUM
#------------------------------------------------
		pdbSum_file_name = pdb_id + "-binding.dat"
		pdbSum_file = open(pdbSum_file_name,"w")
		lowpdb_id= pdb_id.lower() # PDBsum requires lowercase of pdbid
		print lowpdb_id
		url_pdbsum = 'http://www.ebi.ac.uk/thornton-srv/databases/pdbsum/%s/grow.out' % lowpdb_id
		print url_pdbsum
		pdbBind = urllib.urlopen(url_pdbsum).read()
		pdbSum_file.write(pdbBind)
		pdbSum_file.close()

		# parse this file
		monomerpdbsum=[]
		print chainName
		readpdbSum_file = open(pdbSum_file_name, "rU")
		for line in readpdbSum_file:
		    a=line.split()
		    if(a[3]==chainName):
			 monomerpdbsum.append(a[4])
		monomerpdbsumList = unique(monomerpdbsum)
		
	
		residList=[]
		i=0
		for n in range(len(pp_prot)-1):
			residList.append (pp_prot[i].get_id()[1])
			i+=1
			
		print len(monomerpdbsumList)

#		PdbSumList=[]
#		i=0
#		for n in range(len(pp_prot)-1):
#		    print residList[i]
#		    if (residList[i] == monomerpdbsumList):
#			PdbSumList.append('binding')
#			print 	PdbSumList		
#		  	i+=1
		
		return monomerpdbsumList



def msv3d(id):
#------------------------------------------------
## Get data from MSV3D
#------------------------------------------------
		id_file_name = id + "-msv3d.dat"
		id_file = open(id_file_name,"w")
		url_msv3d = 'http://decrypthon.igbmc.fr/msv3d/cgi-bin/humsavar?protein=%s' % id
		print url_msv3d
		msv3d = urllib.urlopen(url_msv3d).read()
		id_file.write(msv3d)
		id_file.close()

