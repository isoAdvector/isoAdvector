#!/bin/bash


NzList=(60)
NxList=(100)
dtList=(0.02)
CoList=(04)

#Domain  dimensions
export L=5
export H=3
export ny=1
export y1=-.5
export y2=.5
#Simulation time
export T=4
#Velocity components
export Ux=1
export Uy=0
export Uz=0.5
#Case location
series=$PWD/$1
#Solver
export solver=scalarTransportFoam


mkdir $series

for n in ${!NzList[*]} 
do
    export nx=${NxList[$n]}
    export nz=${NzList[$n]}
    export dt=${dtList[$n]}
    caseName=N${NzList[$n]}Co${CoList[$n]}
    caseDir=$series/$caseName
    echo $caseDir
    cp -r baseCase $caseDir
    envsubst < baseCase/constant/polyMesh/blockMeshDict > ${caseDir}/constant/polyMesh/blockMeshDict
    envsubst < baseCase/system/controlDict > ${caseDir}/system/controlDict
    envsubst < baseCase/0.org/U > ${caseDir}/0.org/U
    touch ${caseDir}/${caseName}.foam
#    if [ "$solver"="scalarTransportFoam" ];
#    then
#        echo "Turning alpha1 into T since solver is scalarTransportFoam"
#        sed -i 's/\[0 0 0 0 0 0 0\];/\[0 0 0 1 0 0 0\];/' $caseDir/0.org/alpha1
#        sed -i 's/alpha/T/' $caseDir/0.org/alpha1
#        mv $caseDir/0.org/alpha1 $caseDir/0.org/T
#    fi
done
