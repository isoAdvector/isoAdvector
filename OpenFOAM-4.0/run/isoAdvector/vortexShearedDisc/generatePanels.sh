#!/bin/bash

appList=(mulesFoam)
#appList=(isoAdvector)
schemeList=(MULES)
#schemeList=(isoAdvector)
meshList=(hex tri poly)
#meshList=(poly)
CoList=(0.1)
#CoList=(0.5)

pvsmFile=spiralDisc.pvsm
#NzList=(100 200 400)
NzList=(200)

for nn in ${!meshList[*]}
do
    meshType=${meshList[$nn]}

    for mm in ${!appList[*]}
    do
        application=${appList[$mm]}
        scheme=${schemeList[$mm]}

        figDir=figs
        mkdir --parents ./$figDir

        #Case location
        series=$PWD/$scheme/$meshType

        for n in ${!NzList[*]} 
        do
            nz=${NzList[$n]}

            for m in ${!CoList[*]}
            do
                Co=${CoList[$m]}
                caseName=N${nz}Co${Co}
                caseDir=$series/$caseName
                echo $caseDir

                pvbatch --use-offscreen-rendering $FOAM_USER_APPBIN/pvbatchprocess.py \
                        --casedir $scheme/$meshType/$caseName \
                        --statefile $pvsmFile \
                        --outputdir ../../../$figDir \
                        --imagebase ${scheme}_${meshType}_$caseName \
                        --firstTime
            done
        done
    done
done
