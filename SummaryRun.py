#!/usr/bin/env python

########################################
# Written by Eric Thomas & Tyler Kuntz & Nevin Gerek 
# May, 2013
########################################

Usage = ' '

import sys
import SummaryFunctions


if __name__=="__main__" and len(sys.argv)==1:
	print Usage
	sys.exit()

def GetInputFileName(argv, NArg):
  if len(argv) > 2 + NArg:
    print "Please check your command line!"
  else:
    return sys.argv[1]


if __name__=="__main__":

	npBLAST = SummaryFunctions.npBLAST
	unBLAST = SummaryFunctions.unBLAST
	seqBLAST = SummaryFunctions.seqBLAST
	parseBlastFile = SummaryFunctions.parseBlastFile
	matchPDB = SummaryFunctions.matchPDB
	runPDB = SummaryFunctions.runPDB
	runPDBList = SummaryFunctions.runPDBList
        matchPDBHighThroughput = SummaryFunctions. matchPDBHighThroughput
        runPDBListforBA = SummaryFunctions.runPDBListforBA

	# Request user input: choose search by np id or sequence
#	output_name = raw_input("Please enter the name of the output file: \n")	
#	minQueryCoverage = float(raw_input("Please enter the minimum query coverage (in %): \n"))
#	minSeqIdentity = float(raw_input("Please enter the minimum sequence identity (in %): \n"))

# read input file
	output_name = GetInputFileName(sys.argv, 0)
	
	for line in open(output_name ,'rU').readlines():
#	  print line
	  a = line.split()
          minQueryCoverage = float(a[0])
	  minSeqIdentity=float(a[1])


	search_type = float(raw_input("Please select an option:\n(1) Search by NP ID\n(2) Search by Uniprot ID\n(3) Search by sequence \n(4) Put PDB id list \n(5) Blast sequences to match PDBs using NP or Uniprot IDs\no(6) Check pdb for biological assembly \n"))
	while search_type != 1 and search_type != 2 and search_type != 3 and search_type != 4 and search_type != 5 and search_type != 6:
		search_type = float(raw_input("Not a valid selection. Please select an option:\n(1) Search by NP ID \n(2) Search by Uniprot ID \n(3) Search by sequence \n(4) Put PDB id list \n"))
		
	# BLAST search by np id
	if search_type == 1:
		id_file = raw_input("Please enter the name of a .dat file containing NP IDs: \n")
		idList = open(id_file + ".dat", "rU")
		np_id_list = []
		for line in idList:
			np_id = line.split()[0]
			np_id_list.append(np_id)
			print(np_id)
			npBLAST(np_id)
		
		for np_id in np_id_list:
			parseBlastFile(np_id)
			match=matchPDB(np_id, minQueryCoverage, minSeqIdentity)
		
			
        # BLAST search by Uniprot id
        if search_type == 2:
		id_file = raw_input("Please enter the name of a .dat file containing Uniprot IDs: \n")
		idList = open(id_file + ".dat", "rU")
		uniprot_id_list = []
		for line in idList:
			uniprot_id = line.split()[0]
			uniprot_id_list.append(uniprot_id)
			print(uniprot_id)
#			npBLAST(uniprot_id)
			unBLAST(uniprot_id)
			seqBLAST(uniprot_id)
		

		for uniprot_id in uniprot_id_list:
			parseBlastFile(uniprot_id)
			match=matchPDB(uniprot_id, minQueryCoverage, minSeqIdentity)
		

	# BLAST search by sequence
	if search_type == 3:
		single_np_id = raw_input("Please enter the name of a .fasta file containing sequence (name must be NP ID of the sequence): \n")
		np_id_list = [single_np_id]
		seqBLAST(single_np_id)

	# read PDB ids to calculate dfi 
	if search_type == 4:
		pdb_id_file = raw_input("Please enter the name of a .dat file containing PDB IDs: \n")
		pdbidList = open(pdb_id_file + ".dat", "rU")
		pdb_id_list = []
		chain_id_list=[]
		np_id =0
		for line in pdbidList:
			pdb_id = line.split()[0]
			chain_id = line.split()[1]
			pdb_id_list.append(pdb_id)
			chain_id_list.append(chain_id)
			print pdb_id,chain_id
			runPDBList(np_id,pdb_id,chain_id)

	#find matching pubs with blast without calculating StrD parms
	if search_type == 5:
		id_file = raw_input("Please enter the name of a .dat file containing NP or UniProt IDs: \n")
		idList = open(id_file + ".dat", "rU")
		id_list = []
		for line in idList:
			id = line.split()[0]
			id_list.append(id)
			print(id)
			npBLAST(id)
		
		pdb_summary = open(output_name + ".out", "w")
		print >> pdb_summary,"NP \tPDB \t#C \tMR \tQC \tSI"

		for id in id_list:
			parseBlastFile(id)
			match=matchPDBHighThroughput(id, minQueryCoverage, minSeqIdentity)
			if(match is not None):
				pdb_summary.write(match)
			else:
			       print>>pdb_summary,id
		pdb_summary.close()
		
	# read PDB ids if they have BA or not
	if search_type == 6:
		pdb_id_file = raw_input("Please enter the name of a .dat file containing PDB IDs: \n")
		pdbidList = open(pdb_id_file + ".dat", "rU")
		pdb_id_list = []
		chain_id_list=[]
		np_id =0
		for line in pdbidList:
			pdb_id = line.split()[0]
			chain_id = line.split()[1]
			pdb_id_list.append(pdb_id)
			chain_id_list.append(chain_id)
			print pdb_id,chain_id
			runPDBListforBA(np_id,pdb_id,chain_id)
		
