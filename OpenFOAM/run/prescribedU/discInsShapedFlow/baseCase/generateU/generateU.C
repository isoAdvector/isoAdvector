/*---------------------------------------------------------------------------*\
              Original work | Copyright (C) 2016-2017 DHI
              Modified work | Copyright (C) 2016-2017 OpenCFD Ltd.
              Modified work | Copyright (C) 2017-2018 Johan Roenby
-------------------------------------------------------------------------------

License
    This file is part of IsoAdvector, which is an unofficial extension to
    OpenFOAM.

    IsoAdvector is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    IsoAdvector is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with IsoAdvector. If not, see <http://www.gnu.org/licenses/>.

Application
    generateU

Description
    Generates velocity field for interface advection case with a disc deformed
    into an S-shape. Taken from:

    Ahn, Hyung Taek, og Mikhail Shashkov.
    “Adaptive moment-of-fluid method”.
    J. Comput. Phys. 228, nr. 8 (maj 2009): 2792–2821.
    doi:10.1016/j.jcp.2008.12.031.

Author
    Johan Roenby, STROMNING, all rights reserved.

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

    {
        scalarField X = mesh.C().component(0);
        scalarField Z = mesh.C().component(2);
        scalarField u = 0.25*(4.0*X-2.0 + pow(4.0*Z-2.0,3));
        scalarField w = -0.25*(4.0*Z-2.0 + pow(4.0*X-2.0,3));

        //Setting U in all cells
        forAll(U,ci)
        {
            U[ci] = u[ci]*vector(1.0,0.0,0.0) + 0.0*vector(0.0,1.0,0.0) + w[ci]*vector(0.0,0.0,1.0);
        }

    }

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

    {
        scalarField Xf = mesh.Cf().component(0);
        scalarField Zf = mesh.Cf().component(2);
        scalarField uf = 0.25*(4.0*Xf-2.0 + pow(4.0*Zf-2.0,3));
        scalarField wf = -0.25*(4.0*Zf-2.0 + pow(4.0*Xf-2.0,3));
        forAll(phi,fi)
        {
            vector Uf = uf[fi]*vector(1.0,0.0,0.0) + wf[fi]*vector(0.0,0.0,1.0);
            phi[fi] = (Uf & (mesh.Sf()[fi]));
        }
    }

    scalarField sumPhi = fvc::surfaceIntegrate(phi);
    scalar maxMagSumPhi(0.0);
    label maxLabel(0);
    forAll(sumPhi,ci)
    {
        scalar msp = mag(sumPhi[ci]);
        if (msp > maxMagSumPhi )
        {
            maxMagSumPhi = msp;
            maxLabel = ci;
        }
    }
    Info << "maxMagSumPhi/Vi = " << maxMagSumPhi << " at cell " << maxLabel << endl;

    ISstream::defaultPrecision(18);

    U.write();
    phi.write();

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
