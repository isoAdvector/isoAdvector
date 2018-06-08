#!/bin/bash

#Use gmsh to generate a unit square domain tessellated by triangular prisms.
#
#Takes an integer argument which is the number of cells in each direction.
#
#Tested with gmsh 3.0.6 binary downloadedd from http://gmsh.info/#Download
#
#and placed in $FOAM_USER_APPBIN
#
#Johan Roenby, STROMNING, 2018


caseName=triUnitSquare
rm -rf $caseName
cp -r $FOAM_TUTORIALS/basic/potentialFoam/pitzDaily $caseName

echo "nx=${1};" >> $caseName/triSquare.geo

echo 'Point(1) = {-0, -0.005, 0, 1e+22};
Point(2) = {1, -0.005, 0, 1e+22};
Point(3) = {1, -0.005, 1, 1e+22};
Point(4) = {0, -0.005, 1, 1e+22};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(6) = {4, 1, 2, 3};

Plane Surface(6) = {6};
Physical Volume("internal") = {1};

Transfinite Line{1} = nx;
Transfinite Line{2} = nx;
Transfinite Line{3} = nx;
Transfinite Line{4} = nx;

Extrude {0, 0.1, 0} {
  Surface{6};
  Layers{1};
  Recombine;
}

Physical Surface("back") = {6};
Physical Surface("rim") = {27, 15, 19, 23};
Physical Surface("front") = {28};' >> $caseName/triSquare.geo


cd $caseName
touch ${caseName}.foam
gmsh -3 triSquare.geo  > gmsh.log 2>&1
gmshToFoam triSquare.msh > gmshToFoam.log 2>&1

echo "Mesh generated in $caseName/constant/polyMesh"
