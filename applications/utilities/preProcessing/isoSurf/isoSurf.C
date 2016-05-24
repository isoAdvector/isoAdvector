/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of the IsoAdvector source code library, which is an 
	unofficial extension to OpenFOAM.

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
    isoSurf

Description
    Uses isoCutter to create a volume fraction field from either a cylinder, 
    a sphere or a plane.

Author
    Johan Roenby, DHI, all rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "isoCutter.H"
#include "isoCutFace.H"
#include "isoCutCell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Reading field alpha1\n" << endl;
    volScalarField alpha1
    (
        IOobject
        (
            "alpha1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

	Info<< "Reading isoSurfDict\n" << endl;

	IOdictionary isoSurfDict
	(
		IOobject
		(
			"isoSurfDict",
			runTime.system(),
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	);

	const word surfType = isoSurfDict.lookup("type");
	const vector centre = isoSurfDict.lookup("centre");
	const vector direction = isoSurfDict.lookupOrDefault<vector>("direction",vector::zero);
	const scalar radius = isoSurfDict.lookupOrDefault<scalar>("radius",0.0);

	const scalarField x = mesh.points().component(0);
    const scalarField y = mesh.points().component(1);
    const scalarField z = mesh.points().component(2);
    scalar f0 = 0.0;
	scalarField f(x.size());

	if ( surfType == "plane" )
	{
		f = -(mesh.points() - centre) & (direction/mag(direction));
		f0 = 0.0;
	}
	else if ( surfType == "sphere" )
	{
		f = -sqrt(pow((x-centre[0]),2) + pow((y-centre[1]),2) + pow((z-centre[2]),2));
		f0 = -radius;
	}
	else if ( surfType == "cylinder" )
	{
		f = -sqrt(pow(mag(mesh.points()-centre),2) - pow(mag((mesh.points()-centre) & direction),2));
		f0 = -radius;
	}
	else
	{
		Info << "Invalid surface type specified" << endl;
		Info << "Aborting..." << endl;
	}

    //Define function on mesh points and isovalue
	
	//Calculating alpha1 volScalarField from f = f0 isosurface

    //Setting internal alpha1 values
    isoCutCell icc(mesh,f);
    scalarField& alphaIn = alpha1.internalField();
        
    forAll(mesh.cells(),ci)
    {
        const label cellStatus = icc.calcSubCell(ci,f0);
        if (cellStatus != 1) //I.e. if cell not entirely above isosurface
        {
            alphaIn[ci] = icc.VolumeOfFluid();
        }
    }
    
    //Setting boundary alpha1 values
    isoCutFace icf(mesh, f);
    forAll(mesh.boundary(), patchI)
    {
        if (mesh.boundary()[patchI].size() > 0)
        {
            fvPatchScalarField& alphaBP = alpha1.boundaryField()[patchI];

            forAll(alphaBP, faceI)
            {
                const label fLabel = faceI + mesh.boundary()[patchI].start();
                const label faceStatus = icf.calcSubFace(fLabel, f0);
                if (faceStatus != 1) //I.e. if face not entirely above isosurface
                {
                    alphaBP[faceI] = mag(icf.subFaceArea());
                }
            }
        }
    }
	
	ISstream::defaultPrecision(18);
	
    alpha1.write(); //Writing volScalarField alpha1

	Info << "sum(alpha*V) = " << gSum(mesh.V()*alphaIn)
	 << ", 1-max(alpha1) = " << 1 - gMax(alphaIn)
	 << "\t min(alpha1) = " << gMin(alphaIn) << endl;
	
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
