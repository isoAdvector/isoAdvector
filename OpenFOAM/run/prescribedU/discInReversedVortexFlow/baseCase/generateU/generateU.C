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
    Generates velocity field for the classical test case with a sphere 
    deformed into a shape with long tongues and back again.

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

    const scalar pi = Foam::constant::mathematical::pi;

    const scalarField x(mesh.C().component(vector::X));
    const scalarField y(mesh.C().component(vector::Y));
    const scalarField z(mesh.C().component(vector::Z));

    vectorField& Uc = U.primitiveFieldRef();
    Uc.replace(vector::X, -sin(2*pi*z)*sqr(sin(pi*x)));
    Uc.replace(vector::Y, scalar(0));
    Uc.replace(vector::Z, sin(2*pi*x)*sqr(sin(pi*z)));

    Info<< "Reading/calculating face flux field phi\n" << endl;

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
        linearInterpolate(U) & mesh.Sf()
    );
    
    // Calculating incompressible flux based on stream function
    const scalarField xp(mesh.points().component(0));
    const scalarField yp(mesh.points().component(1));
    const scalarField zp(mesh.points().component(2));
    const scalarField psi((1.0/pi)*sqr(sin(pi*xp))*sqr(sin(pi*zp)));

    scalarField& phic = phi.primitiveFieldRef();
    forAll(phic, fi)
    {
        phic[fi] = 0;
        const face& f = mesh.faces()[fi];
        const label nPoints = f.size();

        forAll(f, fpi)
        {
            const label p1 = f[fpi];
            const label p2 = f[(fpi + 1) % nPoints];
            phic[fi] += 0.5*(psi[p1] + psi[p2])*(yp[p2] - yp[p1]);
        }
    }

    surfaceScalarField::Boundary& phibf = phi.boundaryFieldRef();
    forAll(phibf, patchi)
    {
        scalarField& phif = phibf[patchi];
        const label start = mesh.boundaryMesh()[patchi].start();

        forAll(phif, fi)
        {
            phif[fi] = 0;
            const face& f = mesh.faces()[start + fi];
            const label nPoints = f.size();

            forAll(f, fpi)
            {
                const label p1 = f[fpi];
                const label p2 = f[(fpi + 1) % nPoints];
                phif[fi] += 0.5*(psi[p1] + psi[p2])*(yp[p2] - yp[p1]);
            }
        }
    }

    ISstream::defaultPrecision(18);

    U.write();
    phi.write();
    
    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
