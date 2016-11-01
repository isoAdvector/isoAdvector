#!/bin/bash

#Input

#appList=(isoAdvect interFoam passiveAdvectionFoam passiveAdvectionFoam)
appList=(isoAdvector)
#schemeList=(isoAdvector MULES HRIC CICSAM)
schemeList=(isoAdvector)
meshList=(hex tri poly)
CoList=(0.5)

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

#Location of tri meshes
triMeshDir=$ISOADVECTOR_ROOT_DIR/meshes/triangularPrisms_5x01x1

if [ ! -d "$triMeshDir" ];
then
    echo " "
    echo "Warning: "
    echo "The triangular prism meshes cannot be found at the expected location:"
    echo " "
    echo "\$ISOADVECTOR_ROOT_DIR/meshes/triangularPrisms_5x01x1"
    echo " "
    echo "Therefore only hex mesh cases will run properly."
    echo "Please make sure that you have download the meshes."
    echo " "
fi

#Looping through all mesh, solver, and CFL number combinations and generating cases
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
            NzList=(30 60 120)
            NxList=(50 100 200)
        elif [ -d "$triMeshDir" ];
        then
            NzList=(coarse mid fine)
        else
            NzList=()
        fi

        mkdir --parents $series

        for n in ${!NzList[*]}
        do
            for m in ${!CoList[*]}
            do
                Co=${CoList[$m]}
                if [ "$meshType" = "hex" ];
                then
                    caseName=N${NzList[$n]}Co${Co}
                else
                    caseName=${NzList[$n]}Co${Co}
                fi
                caseDir=$series/$caseName
                echo $caseDir
                cp -r baseCase $caseDir

                foamParmSet maxAlphaCo "$Co" $caseDir/system/controlDict
                foamParmSet application "$application" $caseDir/system/controlDict

                if [ "$scheme" = "CICSAM" ];
                then
                    foamParmSet 'div(phi,alpha1)'  "Gauss $scheme 0.5" $caseDir/system/fvSchemes
                fi
                if [ "$scheme" = "HRIC" ];
                then
                    foamParmSet 'div(phi,alpha1)'  "Gauss $scheme" $caseDir/system/fvSchemes
                fi

                if [ "$meshType" = "hex" ];
                then
                    nx=${NxList[$n]}
                    nz=${NzList[$n]}
                    foamParmSet L "$L" $caseDir/constant/polyMesh/blockMeshDict
                    foamParmSet H "$H" $caseDir/constant/polyMesh/blockMeshDict
                    foamParmSet nx "$nx" $caseDir/constant/polyMesh/blockMeshDict
                    foamParmSet nz "$nz" $caseDir/constant/polyMesh/blockMeshDict
                    foamParmSet 'internalField' "uniform (1 0 .5)" $caseDir/0.org/U
                    cp -r $caseDir/0.org $caseDir/0
                else
                    cp $triMeshDir/polyMesh_${NzList[$n]}/* $caseDir/constant/polyMesh/
                    foamParmSet 'internalField' "uniform (1 0 0)" $caseDir/0.org/U
                    cp -r $caseDir/0.org $caseDir/0
                    if [ "$meshType" = "poly" ];
                    then
                        mkdir $caseDir/logs
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

echo " "
echo "Cases generated. To run cases type"
echo " "
echo "foamRunCasesIn " $schemeList

fi
