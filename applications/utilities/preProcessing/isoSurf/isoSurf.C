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
    isoSurf

Description
    Uses isoCutter to create a volume fraction field from either a cylinder, 
    a sphere or a plane.

Author
    Johan Roenby, DHI, all rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "isoCutter.H"

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
    Foam::isoCutter cutter(mesh,f);
    cutter.subCellFractions(f0,alpha1);
	alpha1.correctBoundaryConditions();
	
	ISstream::defaultPrecision(18);
	
    alpha1.write(); //Writing volScalarField alpha1

	Info << "sum(alpha*V) = " << sum(mesh.V()*alpha1).value() 
	 << ", max(alpha1)-1 = " << max(alpha1).value()-1.0
	 << "\t min(alpha1) = " << min(alpha1).value() << endl;
	
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
