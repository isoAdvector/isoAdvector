#!/bin/bash

meshType=tet
dxList=(01)
dtList=(0.01)
CoList=(05)

export T=4
export Ux=0
export Uy=0
export Uz=1
export solver=isoAdvector

for n in ${!dxList[*]} 
do
    export dx=${dxList[$n]}
    export dt=${dtList[$n]}
    caseName=dx${dxList[$n]}Co${CoList[$n]}
    echo $caseName
    cp -r baseCase $caseName
    cp -r ../meshes/$meshType/${dxList[$n]}/* ${caseName}/constant/polyMesh/
    envsubst < baseCase/system/controlDict > ${caseName}/system/controlDict
    envsubst < baseCase/0.org/U > ${caseName}/0.org/U
    touch ${caseName}/${caseName}.foam
done
