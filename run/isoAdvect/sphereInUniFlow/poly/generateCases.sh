#!/bin/bash

meshType=tet
dxList=(01 005 0025)
dtList=(0.01 0.005 0.0025)
CoList=(05 05 05)

export T=4
export Ux=0
export Uy=0
export Uz=1
export solver=isoAdvector

mkdir $solver

for n in ${!dxList[*]} 
do
    export dx=${dxList[$n]}
    export dt=${dtList[$n]}
    caseName=dx${dxList[$n]}Co${CoList[$n]}
    echo $caseName
    caseDir=$solver/$caseName
    cp -r baseCase $caseDir
    cp -r ../meshes/$meshType/${dxList[$n]}/* ${caseDir}/constant/polyMesh/
    envsubst < baseCase/system/controlDict > ${caseDir}/system/controlDict
    envsubst < baseCase/0.org/U > ${caseDir}/0.org/U
    touch ${caseDir}/${caseName}.foam
done
