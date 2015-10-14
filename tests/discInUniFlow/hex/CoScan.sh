#!/bin/bash

#Do "./CoScan generate" to generate cases. Do "./CoScan run" to run cases

#Name of directory (in current directory in which runs are stored)
testName=CoScan_Nx40

#Solver
SOLVERS=(isoAdvect interFoam scalarTransportFoam scalarTransportFoam)
VOFSCHEMES=(none none CICSAM HRIC)
VOFPARMS=(0 0 0.5 0)
NAMES=(isoAdvect MULES CICSAM HRIC)

export nz=120
export nx=200

#dtList=(0.001 0.005 0.025)
dtList=(0.0005 0.0025 0.0125)
CoList=(002 01 05)
writeIntervalList=(25 5 1)

if [ "$1" = "generate" ];
then
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

mkdir $testName

#Making directory for each VOF method
    for m in ${!NAMES[*]}
    do
        mkdir ${testName}/${NAMES[$m]}
    done
fi

#Constructing cases
for n in ${!dtList[*]}
do
    export dt=${dtList[$n]}
	export wInt=${writeIntervalList[$n]}

    for m in ${!SOLVERS[*]}
    do
        export solver=${SOLVERS[$m]}
        export VOFSCHEME=${VOFSCHEMES[$m]}
        export VOFPARM=${VOFPARMS[$m]}
        name=${NAMES[$m]}
        caseName=Co${CoList[$n]}
        caseDir=${testName}/${name}/$caseName
        if [ "$1" = "generate" ];
        then
            echo "Generating case " $caseDir
            cp -r baseCase $caseDir
            envsubst < baseCase/constant/polyMesh/blockMeshDict > ${caseDir}/constant/polyMesh/blockMeshDict
            envsubst < baseCase/system/controlDict > ${caseDir}/system/controlDict
            envsubst < baseCase/system/fvSchemes > ${caseDir}/system/fvSchemes
            envsubst < baseCase/0.org/U > ${caseDir}/0.org/U
            touch ${caseDir}/case.foam
        elif [ "$1" = "run" ];
        then
            echo "Running case " $caseDir
            cd $caseDir
            ./Allrun &
            cd -
        fi
    done
done
