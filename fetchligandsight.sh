#!/bin/bash
#reads in a list or pdbids
#and fetches the ligand binding
#sights. 
#Avishek Kumar avishek.kumar@asu.edu 

#input check for file

PROG=$(basename $0)
LIGFIL="ligands.txt"

usage="
==========================
Get Ligand Binding Sights 
==========================

Usage
-----
${PROG} pdblist.txt 

Description 
-----------
Fetches ligand binding sights and dumps the output into a text file. The
sites are then formatted to be read into dfi in the file $LIGFIL

Example 
-------
pdbid.txt: 
1a2x
1dc2 
1ake 

"

if [ -f $LIGFIL ]
then
    echo "Deleting the $LIGFIL to start anew."
    rm -v $LIGFIL
fi

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
    a=$(awk '!/#/ {print $5 " " $14}' ${f}.out | awk '{print $2$1}')
    echo ${f}:${a} >> ${LIGFIL}
done
