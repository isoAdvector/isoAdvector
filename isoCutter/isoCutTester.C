/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is not part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    isoCutTester

Description
    Testing isoCutter class

Author
    Johan Roenby, DHI

\*---------------------------------------------------------------------------*/

#include "isoCutter.H"
#include "argList.H"
#include "volPointInterpolation.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    //Read in mesh from case directory
    argList::noBanner();
    Foam::argList args(argc, argv);
    Foam::Time runTime(Foam::Time::controlDictName, args);

    Foam::fvMesh mesh
    (
        Foam::IOobject
        (
            Foam::fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );

    volScalarField alpha1
    (
        IOobject
        (
            "alpha1",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    surfaceScalarField alphaf
    (
        IOobject
        (
            "alphaf",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        0*mag(mesh.Sf())/mag(mesh.Sf())
    );

    runTime++;

    //Define function on mesh points and isovalue
    const scalarField x = mesh.points().component(0);
    const scalarField y = mesh.points().component(1);
    const scalarField z = mesh.points().component(2);
    scalar pi = constant::mathematical::pi;
    const scalarField f = .5 + 0.2*sin(2*pi*x)*exp(-pow((y-.5)/.5,2)) - z;
//    const scalarField f = -.26 - .23423*z+.8764*y-.1203*x;

    const scalar f0(0);

	//Calculating alpha1 volScalarField from f = f0 isosurface
    Foam::isoCutter cutter(mesh,f,f0);
    cutter.subCellFractions(f0,alpha1);
	alpha1.write(); //Writing volScalarField alpha1
    cutter.write(); //Writing cutCells to ply files

	//Interpolating alpha1 to vertices
	volPointInterpolation vpi(mesh);
	alpha1.correctBoundaryConditions();
	scalarField alphap = vpi.interpolate(alpha1);

	//Finding for each cell the isoValue that cuts the cell to the right volume fraction
    Foam::isoCutter cutter2(mesh,alphap,f0);
	volScalarField isoVals = 0*alpha1;
	scalar tol = 1e-6;
	label maxIter = 500;
	cutter2.vofCutCells(alpha1, tol, maxIter, isoVals);
//    cutter.subFaceFractions(f0,alphaf);
//    alphaf.write();


    return 0;
}
// ************************************************************************* //