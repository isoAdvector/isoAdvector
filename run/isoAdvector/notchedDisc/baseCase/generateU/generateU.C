/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    generateU

Description
    Generates velocity field for the classical test case with a notched disc
    in a rigid body rotation flow.

Author
    Johan Roenby, DHI, all rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Reading field U\n" << endl;

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh
    );

    Info<< "Calculating field U\n" << endl;
    {
        const volVectorField& C = mesh.C();
        forAll(U, ci)
        {
            scalar X = C[ci].component(vector::X);
            scalar Z = C[ci].component(vector::Z);
            scalar u = -2*M_PI*(Z-.5);
            scalar w = 2*M_PI*(X-.5);
            U[ci] =  u*vector(1,0,0) + w*vector(0,0,1);
        }
    }

    Info<< "Reading field phi\n" << endl;

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        0*(linearInterpolate(U) & mesh.Sf())
    );

    Info<< "Calculating field phi\n" << endl;

    const surfaceVectorField& Cf = mesh.Cf();
    const surfaceVectorField& Sf = mesh.Sf();
    forAll(phi, fi)
    {
        scalar X = Cf[fi].component(vector::X);
        scalar Z = Cf[fi].component(vector::Z);
        scalar u = -2*M_PI*(Z-.5);
        scalar w = 2*M_PI*(X-.5);
        vector Uf =  u*vector(1,0,0) + w*vector(0,0,1);
        phi[fi] = Uf & Sf[fi];
    }

    surfaceScalarField::GeometricBoundaryField& phip = phi.boundaryField();

    forAll(phip, patchi)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchi];
        if
        (
            !isA<processorPolyPatch>(pp)
         && !isA<emptyPolyPatch>(pp)
        )
        {
            fvsPatchScalarField& phib = phip[patchi];
            const fvsPatchVectorField& Cf = mesh.Cf().boundaryField()[patchi];
            const fvsPatchVectorField& Sf = mesh.Sf().boundaryField()[patchi];

            forAll(phib, fi)
            {
                scalar X = Cf[fi].component(vector::X);
                scalar Z = Cf[fi].component(vector::Z);
                scalar u = -2*M_PI*(Z-.5);
                scalar w = 2*M_PI*(X-.5);
                vector Uf =  u*vector(1,0,0) + w*vector(0,0,1);
                phib[fi] = Uf & Sf[fi];
            }
        }
    }

/*
    //Checking that phi's of a cell sum to zero
    //For the polygonal meshes there are some strange continuity errors but
    //these only seem to appear at the boundary and are therefore irrelevant
    //for this test case.
    volScalarField sumPhiByV
    (
        IOobject
        (
            "sumPhiByV",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::surfaceIntegrate(phi)
    );

    sumPhiByV = fvc::surfaceIntegrate(phi);
    Info << "max(sum(phi)/V) = " << max(sumPhiByV) << endl;
    sumPhiByV.write();
*/

    Info<< "Writing U and phi\n" << endl;
    phi.write();
    U.write();

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
