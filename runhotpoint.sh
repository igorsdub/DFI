#!/bin/bash
usage="
=============
Run Hotpoints 
=============

Description
-----------
Runs batch job for hotpoints. Need pdbs with chains AB

Usage
-----
./runhotpoints.sh \$(cat list.txt)
./runhotpoint.sh 1a2x.pdb 

Example
--------

list.txt:
1a2x.pdb 
1efv.pdb 

Output 
------
batchlist.txt 

"

#input check 
if [ -z $1 ]
then
    echo "ERROR:No input"
    echo "$usage"
    exit 1 
fi 


if [ -f batchlist.dat ]
then
    echo "Removing old batchlist.txt file ..."
    rm -vf batchlist.dat 
fi


for pdb in $*; 
do
    echo $pdb
   
    PDBF=${pdb}
    if !([ -f ${PDBF} ])
    then
	echo "${PDBF} is not a file downloading from pdb."
	./download_pdb.sh ${PDBF}
	if !([ -f ${PDBF} ])
	then
	    echo "${PDBF} not in PDB"
	    echo "usage"
	    exit 1 
	fi
    fi

    interface="${PDBF%.pdb}CD"
    #Runhotpoints for chains AB 
    ./predict_hotspot.py $interface ./ ./ ./ ./
    a=$(grep ",H$" ${interface}.output | awk -F, '{print $3$1}')
    echo "${PDBF}:"${a}"" 
    echo "${PDBF}:"${a}"" >> batchlist.dat 
done
exit 0 

