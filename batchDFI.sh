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
Coming soon: Make individual directories for info.  

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

FILE=${1}
if [ -z ${FILE} ];
then
    die "No input file" 2 
fi

cat ${FILE} | while read f; 
do
    echo $f; 
    pdbfil=$(echo $f | awk -F: '{print $1}')
    chain=$(echo $f | awk -F: '{print $2}')
    ./dfi.py --pdb $pdbfil --fdfi $chain 
done
