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

#include "timeSelector.H"
#include "fvCFD.H"
#include "isoCutter.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
	timeSelector::addOptions();

    #include "setRootCase.H"
    #include "createTime.H"
    #include "readTimeControls.H"
	instantList instList = runTime.times();
	//Setting time to first time which is not 'constant'
	runTime.setTime(instList[1],0); 	
    #include "createMesh.H"
	
	scalar delta0 = Foam::pow(sum(mesh.V()).value(),1.0/3.0);
//	Info << "Average edge lengt = " << delta0 << endl;
	
	Info<< "Reading field alpha1\n" << endl;
	volScalarField alpha1
	(
		IOobject
		(
			"alpha1",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		),
		mesh
	);
	scalar V0 = sum(mesh.V()*alpha1).value();
	Info << "Initial water volume V0 = " << V0 << endl;
	
	Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

	const vector U0 = U[0];
	Info << "Velocity U0 = " << U0 << endl;

    volScalarField alpha1Exact = alpha1; //Values are overwritten
	
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
	const vector centre0 = isoSurfDict.lookup("centre");
	const vector direction = isoSurfDict.lookupOrDefault<vector>("direction",vector::zero);
	const scalar radius = isoSurfDict.lookupOrDefault<scalar>("radius",0.0);

	const scalarField x = mesh.points().component(0);
    const scalarField y = mesh.points().component(1);
    const scalarField z = mesh.points().component(2);
    scalar f0 = 0.0;
	scalarField f(x.size());

	volPointInterpolation vpi(mesh);
	volScalarField alpha101 = alpha1; //Convenience copy: values are overwritten in next line
	volScalarField alpha199 = alpha1; //Convenience copy: values are overwritten in next line
	volScalarField alpha1Exact01 = alpha1; //Convenience copy: values are overwritten in next line
	volScalarField alpha1Exact99 = alpha1; //Convenience copy: values are overwritten in next line

	//Time loop
    instantList timeDirs = timeSelector::select0(runTime, args);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
//        #include "readTimeControls.H"
	
//		Info << "Before runTime++: runTime.timeName() = " << runTime.timeName() << endl;
//		Info << "Before runTime++: runTime.time().value() = " << runTime.time().value() << endl;
//        runTime++;
//		Info << "After runTime++: runTime.timeName() = " << runTime.timeName() << endl;
//		Info << "After runTime++: runTime.time().value() = " << runTime.time().value() << endl;

		volScalarField alpha1
		(
			IOobject
			(
				"alpha1",
				runTime.timeName(),
				mesh,
				IOobject::MUST_READ,
				IOobject::NO_WRITE
			),
			mesh
		);
		
		vector centre = centre0 + U0*runTime.time().value();
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
		cutter.subCellFractions(f0,alpha1Exact);
		scalar E1 = sum(mag(alpha1Exact-alpha1)*mesh.V());
		vector cmExact = sum(alpha1Exact*mesh.V()*mesh.C())/sum(mesh.V());
		vector cm = sum(alpha1*mesh.V()*mesh.C())/sum(mesh.V());
		
		scalarField ap = vpi.interpolate(alpha1);
		Foam::isoCutter cutter2(mesh,ap);
		
		//Volume if alpha = 0.1 isosurface for calculated solution
		cutter2.subCellFractions(.01,alpha101);	
		scalar V01 = sum(mesh.V()*alpha101).value();

		//Volume if alpha = .99 isosurface for calculated solution		
		cutter2.subCellFractions(.99,alpha199);	
		scalar V99 = sum(mesh.V()*alpha199).value();

		scalarField apExact = vpi.interpolate(alpha1Exact);
		Foam::isoCutter cutter3(mesh,apExact);
		
		//Volume if alpha = 0.1 isosurface for exact solution
		cutter3.subCellFractions(.01,alpha1Exact01);	
		scalar V01Exact = sum(mesh.V()*alpha1Exact01).value();

		//Volume if alpha = .99 isosurface for exact solution		
		cutter3.subCellFractions(.99,alpha1Exact99);	
		scalar V99Exact = sum(mesh.V()*alpha1Exact99).value();
		
//		scalar PI = Foam::constant::mathematical::pi;
//		scalar delta = Foam::pow(3.0/(4.0*PI)*V01,1.0/3.0) - Foam::pow(3.0/(4.0*PI)*V99,1.0/3.0);
//		scalar deltaExact = Foam::pow(3.0/(4.0*PI)*V01Exact,1.0/3.0) - Foam::pow(3.0/(4.0*PI)*V99Exact,1.0/3.0);
		scalar dVExact = V01Exact - V99Exact;
		scalar dV = V01 - V99;
		ISstream::defaultPrecision(18);

		Info << "Time: " << runTime.time().value() << endl;
		Info << "E1 = " << E1 << "m3, E1/V0 = " << E1/V0 << endl;
		Info << "(V-V0)/V0 = " << (sum(mesh.V()*alpha1).value()-V0)/V0 << endl;
		Info << "min(alpha1) = " << min(alpha1).value()
			 << ",\t1-max(alpha1) = " << 1-max(alpha1).value() << endl;
		Info << "(dV-dVEx)/dVEx = " << (dV-dVExact)/dVExact << ", where dV = " 
			 << dV << " m3 and dVEx = " << dVExact << " m3" << endl;
		Info << "dCM = " << mag(cm-cmExact) << endl;
//			<< ",\td_Ex/d0 = " << deltaExact/delta0 << endl;
//			<< ",\t(V01-Ex)/Ex = " << (V01-V01Exact)/V01Exact
//			<< ",\t(V99-Ex)/Ex = " << (V99-V99Exact)/V99Exact << endl;		 
	}

    return 0;
}


// ************************************************************************* //
