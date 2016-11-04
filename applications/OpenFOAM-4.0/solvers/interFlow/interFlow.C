/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    interFoam

Description
    Solver for 2 incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

    For a two-fluid approach see twoPhaseEulerFoam.
    
Author
    Johan Roenby, DHI, all rights reserved.


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "isoAdvection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createRDeltaT.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "correctPhi.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    scalar executionTime = runTime.elapsedCpuTime();
    scalar advectionTime = 0;

    while (runTime.run())
    {
        #include "readTimeControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }
        
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            //Advance alpha1 from time t to t+dt
            scalar advectionStartTime = runTime.elapsedCpuTime();
            advector.advect();
            advectionTime += (runTime.elapsedCpuTime() - advectionStartTime);
//            #include "alphaControls.H"
//            #include "alphaEqnSubCycle.H"

            //Clip and snap alpha1 to ensure strict boundedness to machine precision
            alpha1Org = alpha1;
            if ( clipAlphaTol > 0.0 )
            {
                alpha1 = alpha1*
                    pos(alpha1-clipAlphaTol)*neg(alpha1-(1.0-clipAlphaTol)) 
                    + pos(alpha1-(1.0-clipAlphaTol));
            }
            if ( snapAlpha )
            {
                alpha1 = min(1.0,max(0.0,alpha1));
            }

            rho == alpha1*rho1 + (scalar(1) - alpha1)*rho2;
            rhoPhi = advector.getRhoPhi(rho1, rho2);
            
            mixture.correct();

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        alpha1 = alpha1Org;
        const scalar V = gSum(mesh.V().field()*alpha1.internalField());
        Info << "t = " << runTime.time().value() << ",\t sum(alpha*V) = " << V
             << ",\t dev = " << 100*(1.0 - V/V0) << "%"
             << ",\t 1-max(alpha1) = " << 1 - gMax(alpha1.internalField())
             << ",\t min(alpha1) = " << gMin(alpha1.internalField()) 
             << endl;

        if (printSurfCells)
        {
            advector.getSurfaceCells(surfCells);
        }

       if (printBoundCells)
        {
            advector.getBoundedCells(boundCells);
        }

        contErr = fvc::surfaceIntegrate(phi)*runTime.deltaT();

        runTime.write();

        scalar newExecutionTime = runTime.elapsedCpuTime();
        Info<< "ExecutionTime = " << newExecutionTime << " s,"
            << " ClockTime = " << runTime.elapsedClockTime() << " s,"
            << " timeStepTime = " << newExecutionTime - executionTime << " s,"
            << " advection fraction: " << 100*advectionTime/newExecutionTime
            << "%" << nl << endl;
        executionTime = runTime.elapsedCpuTime();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //