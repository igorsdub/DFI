#!/bin/bash

pgrmname=$(basename $0)

usage="
batchUniproPDB.sh
==================

Description
------------
Runs Batch run for blasting UniproID(s) and getting PDBid

Usage
-----
$pgrmname listfile 

"

function die() {
local errmsg="$1" errcode="${2:-1}"
echo "ERROR: ${errmsg}"
echo "${usage}"
exit ${errcode}
}

listfile=${1}
if [ -z $listfile ];
then
    die "No input list" 2
fi

while read f; 
do 
    ./UniproBlastToPdb.py $f; 
done < $listfile 
