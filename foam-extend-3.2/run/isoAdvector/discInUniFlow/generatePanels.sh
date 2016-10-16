#!/bin/bash

appList=(isoAdvector interFoam passiveAdvectionFoam passiveAdvectionFoam)
#appList=(isoAdvector)
schemeList=(isoAdvector MULES HRIC CICSAM)
#schemeList=(isoAdvector)
meshList=(hex tri poly)
#meshList=(hex)
CoList=(0.1 0.2 0.5)
#CoList=(0.1)


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

        if [ "$meshType" = "hex" ];
        then
            NzList=(30 60 120)
        else
            NzList=(20 40 80)
        fi

        for n in ${!NzList[*]} 
        do
            nz=${NzList[$n]}

            for m in ${!CoList[*]}
                do
                    Co=${CoList[$m]}
                    caseName=N${nz}Co${Co}
                    caseDir=$series/$caseName
                    echo $caseDir

                    if [ "$meshType" = "hex" ];
                    then
                        pvsmFile=panel.pvsm
                    else
                        pvsmFile=panelNonhex.pvsm
                    fi

                    pvbatch --use-offscreen-rendering $FOAM_USER_APPBIN/pvbatchprocess.py \
                            --casedir $scheme/$meshType/$caseName \
                            --statefile $pvsmFile \
                            --outputdir ../../../$figDir \
                            --imagebase ${scheme}_${meshType}_$caseName \
                            --latestTime

			done
		done
	done
done
