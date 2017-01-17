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

		//Setting U on all boundary patches
		forAll(mesh.boundary(), patchI)
		{
			const vectorField& Cf = mesh.boundary()[patchI].Cf();
			scalarField Xb = Cf.component(0);
			scalarField Zb = Cf.component(2);
			scalarField ub = 0.25*(4.0*Xb-2.0 + pow(4.0*Zb-2.0,3));
			scalarField wb = -0.25*(4.0*Zb-2.0 + pow(4.0*Xb-2.0,3));
			vectorField& Ub = U.boundaryFieldRef()[patchI];
			forAll(Ub,fi)
			{
				Ub[fi] = ub[fi]*vector(1.0,0.0,0.0) + 0.0*vector(0.0,1.0,0.0) + wb[fi]*vector(0.0,0.0,1.0);
			}
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
