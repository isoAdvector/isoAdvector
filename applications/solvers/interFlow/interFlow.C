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
    interFlow

Description
    Solver for 2 incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

    For a two-fluid approach see twoPhaseEulerFoam.

    This solver is essentially the interFoam solver with MULES replaced by
    IsoAdvector for the interface advection step.

Author
    Johan Roenby, DHI, all rights reserved

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "MULES.H"
#include "subCycle.H"
#include "interfaceProperties.H"
#include "twoPhaseMixture.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"
#include "isoAdvection.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "readTimeControls.H"
    #include "correctPhi.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    scalar executionTime = runTime.elapsedCpuTime();
    scalar advectionTime = 0;
    
    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        scalar advectionStartTime = runTime.elapsedCpuTime();

        twoPhaseProperties.correct();

        //Advance alpha1 from time t to t+dt
        advector.advect();
        
        if (printSurfCells)
        {
            advector.getSurfaceCells(surfCells);
        }
        if (printBoundCells)
        {
            advector.getBoundedCells(boundCells);
        }
        
        Info << "1-max(alpha1) = " << 1-gMax(alpha1.internalField()) 
            << " and min(alpha1) = " << gMin(alpha1.internalField()) << endl;

        //Clip and snap alpha1 to ensure strict boundedness to machine precision
        if ( clipAlphaTol > 0.0 )
        {
            alpha1 = alpha1*pos(alpha1-clipAlphaTol)*neg(alpha1-(1.0-clipAlphaTol)) + pos(alpha1-(1.0-clipAlphaTol));
        }
        if ( snapAlpha )
        {
            alpha1 = min(1.0,max(0.0,alpha1));
        }

        rho == alpha1*rho1 + (scalar(1) - alpha1)*rho2;
        rhoPhi = advector.getRhoPhi(rho1, rho2);
        
//        #include "alphaEqnSubCycle.H"
        interface.correct();

        advectionTime += (advectionStartTime - runTime.elapsedCpuTime());
        
        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
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

        contErr = fvc::surfaceIntegrate(phi)*runTime.deltaT();

        runTime.write();

        scalar newExecutionTime = runTime.elapsedCpuTime();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << "  timeStepTime = " << newExecutionTime - executionTime << " s"
            << "  advection fraction: " << 100*advectionTime/newExecutionTime
            << "%" << nl << endl;
        executionTime = runTime.elapsedCpuTime();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
