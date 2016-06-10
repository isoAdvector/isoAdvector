#!/bin/bash
schemeList=(isoAdvector)

for nn in ${!schemeList[*]}
do
    scheme=${schemeList[$nn]}
    for n in $(ls -tr $scheme)
    do
        caseDir=$scheme/$n
        echo $caseDir
        cd $caseDir
        ./Allrun
        cd -
    done
done
