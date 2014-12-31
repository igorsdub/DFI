#!/bin/bash
#reads in a list or pdbids
#and fetches the ligand binding
#sights. 
#Avishek Kumar avishek.kumar@asu.edu 

#input check for file

PROG=$(basename $0)

usage="
==========================
Get Ligand Binding Sights 
==========================

Usage
-----
${PROG} pdblist.txt 

Description 
-----------
Fetches ligand binding sights and dumps the output into a text file. 

Example 
-------
pdbid.txt: 
1a2x
1dc2 
1ake 

"



if [ ! -f "$1" ]
then
    echo "ERROR: FILE NOT ENTERED"
    echo "$usage"
fi

pdbids=${1}

cat $pdbids | while read f; 
do
    echo "PDBID: $f"
    URL="http://www.ebi.ac.uk/thornton-srv/databases/pdbsum/${f}/grow.out"
    curl ${URL} > ${f}.out 
done




