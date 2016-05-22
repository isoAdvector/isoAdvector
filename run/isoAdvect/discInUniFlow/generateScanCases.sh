#!/bin/bash

appList=(isoAdvector interFoam passiveAdvectionFoam passiveAdvectionFoam)
#appList=(isoAdvector)
schemeList=(isoAdvector MULES HRIC CICSAM)
#schemeList=(isoAdvector)
meshList=(hex tri poly)
#meshList=(hex)
CoList=(0.1 0.2 0.5)

#Location of tri meshes
triMeshDir=../triMeshes

for nn in ${!meshList[*]}
do
	meshType=${meshList[$nn]}
	for mm in ${!appList[*]}
	do
		application=${appList[$mm]}
		scheme=${schemeList[$mm]}

		#Case location
		series=$PWD/$scheme/$meshType

		if [ "$meshType" = "hex" ];
		then
			#Domain  dimensions
			L=5
			H=3
			#Vertical velocity component
			Uz=0.5
			NzList=(30 60 120)
			NxList=(50 100 200)
		else
			NzList=(20 40 80)
			Uz=0.0
		fi

		mkdir --parents $series

		for n in ${!NzList[*]} 
		do
			for m in ${!CoList[*]}
			do
				Co=${CoList[$m]}
				caseName=N${NzList[$n]}Co${Co}
				caseDir=$series/$caseName
				echo $caseDir
				cp -r baseCase $caseDir
				
				./ofset maxAlphaCo "$Co" $caseDir/system/controlDict
				./ofset application "$application" $caseDir/system/controlDict
				./ofset Uz "$Uz" $caseDir/0.org/U

				if [ "$scheme" = "CICSAM" ];
				then
					./ofset 'div(phi,alpha1)'  "Gauss $scheme 0.5" $caseDir/system/fvSchemes
				fi
				if [ "$scheme" = "HRIC" ];
				then
					./ofset 'div(phi,alpha1)'  "Gauss $scheme" $caseDir/system/fvSchemes
				fi

				#Generating mesh
				mkdir $caseDir/logs
				cp -r $caseDir/0.org $caseDir/0
				touch $caseDir/case.foam
				if [ "$meshType" = "hex" ];
				then
					nx=${NxList[$n]}
					nz=${NzList[$n]}
					./ofset L "$L" $caseDir/constant/polyMesh/blockMeshDict
					./ofset H "$H" $caseDir/constant/polyMesh/blockMeshDict
					./ofset nx "$nx" $caseDir/constant/polyMesh/blockMeshDict
					./ofset nz "$nz" $caseDir/constant/polyMesh/blockMeshDict
					blockMesh -case $caseDir > $caseDir/logs/blockMesh.log 2>&1
				else
					cp $triMeshDir/N${NzList[$n]}/* $caseDir/constant/polyMesh/
					if [ "$meshType" = "poly" ];
					then
						#Convert from tet to poly mesh
						polyDualMesh -case $caseDir -overwrite 160 > $caseDir/logs/polyDualMesh.log 2>&1
						#Remove backmost  part of cells
						topoSet -case $caseDir > $caseDir/logs/topoSet.log 2>&1
						subsetMesh -case $caseDir -overwrite c0 -patch front > $caseDir/logs/subsetMesh.log 2>&1
						rm -rf *.obj			
					fi
				fi
				
			done
		done
	done
done
