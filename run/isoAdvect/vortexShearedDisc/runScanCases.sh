#!/bin/bash

#schemeList=(isoAdvector HRIC CICSAM MULES)
schemeList=(isoAdvector)
meshList=(hex tri poly)

for nn in ${!schemeList[*]}
do
	scheme=${schemeList[$nn]}

	for mm in ${!meshList[*]}
	do
		mesh=${meshList[$mm]}
		casesDir=$scheme/$mesh
	
		for n in $(ls $casesDir)
		do
			caseDir=$casesDir/$n
			echo $caseDir
			cd $caseDir
			./Allrun
			cd -
		done
	done
done
