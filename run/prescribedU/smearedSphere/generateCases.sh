#!/bin/bash

#Input

appList=(interFlow)
schemeList=(interFlow)
#NxList=(64)
NxList=(64 128 256)
#dtList=(0.004)
dtList=(0.004 0.002 0.001)
Co=0.5

#End of input

#Check if ISOADVECTOR_ROOT_DIR is set
if [ -z "$ISOADVECTOR_ROOT_DIR" ];
then
    echo " "
    echo "Warning: "
    echo "Please set and export ISOADVECTOR_ROOT_DIR to the "
    echo "root directory of the isoAdvector source code."
    echo "Aborting "
    echo " "
else

# Source utilities
. $ISOADVECTOR_ROOT_DIR/bin/dhiFoamTools

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

        foamParmSet "application" "$application" $caseDir/system/controlDict
        foamParmSet "deltaT" "$dt" $caseDir/system/controlDict

        if [ "$scheme" = "CICSAM" ];
        then
            foamParmSet 'interfaceMethod'  "fvSchemes" $caseDir/system/fvSolution
            foamParmSet 'div(phi,alpha.water)'  "Gauss $scheme 0.5" $caseDir/system/fvSchemes
        fi
        if [ "$scheme" = "HRIC" ];
        then
            foamParmSet 'interfaceMethod'  "fvSchemes" $caseDir/system/fvSolution
            foamParmSet 'div(phi,alpha.water)'  "Gauss $scheme" $caseDir/system/fvSchemes
        fi

        if [ "$scheme" = "MULES" ];
        then
            foamParmSet 'interfaceMethod'  "MULES" $caseDir/system/fvSolution
        fi
        
        nx=${NxList[$n]}
        foamParmSet "nx" "$nx" $caseDir/constant/polyMesh/blockMeshDict
        foamParmSet "ny" "$nx" $caseDir/constant/polyMesh/blockMeshDict
        foamParmSet "nz" "$nx" $caseDir/constant/polyMesh/blockMeshDict
    done
done

echo " "
echo "Cases generated. To run cases type"
echo " "
echo "foamRunCasesIn " $schemeList
echo " "

fi
