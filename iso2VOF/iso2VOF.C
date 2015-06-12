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
    iso2Vof

Description
    Calculates VOF field from the cell cuttings of an isosurface.

Author
	Johan Roenby, DHI, all rights reserved.

\*---------------------------------------------------------------------------*/

#include "isoCutter.H"
#include "argList.H"

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

	#include "isoFun.H"

	//Calculating alpha1 volScalarField from f = f0 isosurface
    Foam::isoCutter cutter(mesh,f);
    cutter.subCellFractions(f0,alpha1);
	alpha1.correctBoundaryConditions();
    alpha1.write(); //Writing volScalarField alpha1
//    cutter.write(); //Writing cutCells to ply files

    return 0;
}
// ************************************************************************* //
