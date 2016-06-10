#!/bin/bash

appList=(isoAdvector)
schemeList=(isoAdvector)
#NxList=(64)
NxList=(64 128 256)
#dtList=(0.004)
dtList=(0.004 0.002 0.001)
Co=0.5

for mm in ${!appList[*]}
do
    application=${appList[$mm]}
    scheme=${schemeList[$mm]}

    #Case location
    series=$PWD/$scheme
    mkdir --parents $series

    for n in ${!NxList[*]}
    do
        Nx=${NxList[$n]}
        dt=${dtList[$n]}
        caseName=N${Nx}Co${Co}
        caseDir=$series/$caseName
        echo $caseDir
        cp -r baseCase $caseDir

        ./ofset application "$application" $caseDir/system/controlDict
        ./ofset deltaT "$dt" $caseDir/system/controlDict

        if [ "$scheme" = "CICSAM" ];
        then
            ./ofset 'div(phi,alpha1)'  "Gauss $scheme 0.5" $caseDir/system/fvSchemes
        fi
        if [ "$scheme" = "HRIC" ];
        then
            ./ofset 'div(phi,alpha1)'  "Gauss $scheme" $caseDir/system/fvSchemes
        fi
        nx=${NxList[$n]}
        ./ofset nx "$nx" $caseDir/constant/polyMesh/blockMeshDict
        ./ofset ny "$nx" $caseDir/constant/polyMesh/blockMeshDict
        ./ofset nz "$nx" $caseDir/constant/polyMesh/blockMeshDict
    done
done
