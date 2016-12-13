/*---------------------------------------------------------------------------*\

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
    along with IsoAdvector.  If not, see <http://www.gnu.org/licenses/>.

Application
    isoSurf

Description
    Uses isoCutter to create a volume fraction field from either a cylinder,
    a sphere or a plane.

Author
    Johan Roenby, DHI, all rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
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

	const scalarField x(mesh.points().component(0));
    const scalarField y(mesh.points().component(1));
    const scalarField z(mesh.points().component(2));
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
	else if ( surfType == "sin" )
	{
        const scalar lambda = isoSurfDict.lookupOrDefault<scalar>("lambda",1);
        const scalar amplitude = isoSurfDict.lookupOrDefault<scalar>("amplitude",.1);
        const vector up = isoSurfDict.lookupOrDefault<vector>("up",vector::zero);
        const scalarField xx((mesh.points()-centre) & direction/mag(direction));
        const scalarField zz((mesh.points()-centre) & up/mag(up));
		f = amplitude*Foam::sin(2*M_PI*xx/lambda) - zz;
		f0 = 0;
	}
	else
	{
		Info << "Invalid surface type specified" << endl;
		Info << "Aborting..." << endl;
	}


    Info << "surfType = " << surfType << endl;

    //Define function on mesh points and isovalue

	//Calculating alpha1 volScalarField from f = f0 isosurface
    isoCutCell icc(mesh,f);
    icc.VolumeOfFluid(alpha1, f0);

	ISstream::defaultPrecision(18);

    alpha1.write(); //Writing volScalarField alpha1

    const scalarField& alpha = alpha1.internalField();
	Info << "sum(alpha*V) = " << gSum(mesh.V()*alpha)
	 << ", 1-max(alpha1) = " << 1 - gMax(alpha)
	 << "\t min(alpha1) = " << gMin(alpha) << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
