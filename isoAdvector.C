/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    interFoam

Description
    Solver for 2 incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

    For a two-fluid approach see twoPhaseEulerFoam.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "MULES.H"
#include "subCycle.H"
#include "interfaceProperties.H"
#include "twoPhaseMixture.H"
#include "turbulenceModel.H"
//#include "pimpleControl.H"
#include "fvIOoptionList.H"
#include "volPointInterpolation.H"
#include "isoCutter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

//    pimpleControl pimple(mesh);

    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "readTimeControls.H"
//    #include "correctPhi.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
	
    while (runTime.run())
    {
        #include "readTimeControls.H"
//        #include "CourantNo.H"
//        #include "alphaCourantNo.H"
//        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

		//Alpha loop
		surfaceScalarField dVf = 0*phi*runTime.deltaT()/mesh.V()[0];
		scalarField& dVtest = dVf;
		
		volPointInterpolation vpi(mesh);
		scalarField alphap = vpi.interpolate(alpha1);
		Foam::isoCutter cutter2(mesh,alphap,0.5); //isoValue not used for anything
		cutter2.updateAlpha(alpha1, phi, U, runTime.deltaT().value(), dVtest);
		volScalarField dV = fvc::surfaceIntegrate(dVf); //For each cell sum contributions from faces with pos sign for owner and neg sign for neighbour (as if it is a flux) and divide by cell volume
		
		forAll(alpha1,ci)
		{
			alpha1[ci] += dV[ci];
		}
		alpha1.correctBoundaryConditions();

		Info << "sum(alpha*V) = " << sum(mesh.V()*alpha1).value() 
			 << ", max(alpha1) = " << max(alpha1).value() 
			 << "\t min(alpha1) = " << min(alpha1).value() << endl;
//			 << "\t min(dt) = " << dt.value() << endl; 

        alpha1.write();
        U.write();
		cutter.write();
		
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
