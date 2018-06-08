#!/bin/bash

#Use gmsh to generate a unit cube domain tessellated by tetrahedrons.
#
#Takes an integer argument which is the number of cells in each direction.
#
#Tested with gmsh 3.0.6 binary downloadedd from http://gmsh.info/#Download
#
#and placed in $FOAM_USER_APPBIN
#
#Johan Roenby, STROMNING, 2018

caseName=tetUnitCube
rm -rf $caseName
cp -r $FOAM_TUTORIALS/basic/potentialFoam/pitzDaily $caseName

echo "nx=${1};" >> $caseName/tetBlock.geo

echo 'SetFactory("OpenCASCADE");
//+
Block(1) = {0, 0, 0, 1, 1, 1};
//+
Transfinite Line {11, 9, 12, 10} = nx Using Progression 1;
//+
Transfinite Line {7, 6, 5, 8, 3, 4, 1, 2} = nx Using Progression 1;
//+
Physical Surface("wall") = {6, 2, 4, 3, 5, 1};
Physical Volume("internal") = {1};' >> $caseName/tetBlock.geo


cd $caseName
touch ${caseName}.foam
gmsh -3 tetBlock.geo  > gmsh.log 2>&1
gmshToFoam tetBlock.msh > gmshToFoam.log 2>&1

echo "Mesh generated in $caseName/constant/polyMesh"
