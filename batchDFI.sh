#!/bin/bash 
#Avishek Kumar avishek.kumar@asu.edu 
#Run bash DFI commands 

PROG=$(basename $0)

usage="
===============
Batch DFI Files
===============

Usage
-----
${PROG} batchlist.txt

Description
-----------
Runs multiple dfi files. 
And dumps the results into a directory with the name of the pdb file.   

Example
-------
Example batchlist.txt file:
1a2x_BA_1.pdb:A15 A95 A98 A101 A102 A118 A119 A126 B17 B20 B21 B22 B24 B29
1dc2.pdb:
"


function die() {
local errmsg="$1" errcode="${2:-1}"
echo "ERROR: ${errmsg}"
echo "${usage}"
exit ${errcode}
}

function check_dir() {
"check if the directory exits and create it if it does not"
if [ ! -d $1 ];
then
    echo "${1} directory does not exist..creating it now"
    mkdir -vp ${1}
fi
}

FILE=${1}
if [ -z ${FILE} ];
then
    die "No input file" 2 
fi

cat ${FILE} | while read f; 
do
    echo $f; 
    pdbfil=$(echo $f | awk -F: '{print $1}')
    pdbid=${pdbfil%".pdb"}
    echo $pdbid
    chain=$(echo $f | awk -F: '{print $2}')
    ./dfi.py --pdb $pdbfil --fdfi $chain 
    if [ "$?" -eq "0" ]
    then
	check_dir $pdbid
	mv -v ${pdbid}-* $pdbid
    fi
done
