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
    interFlow

Description
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
#include "isoAdvection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "readGravitationalAcceleration.H"
#   include "readPIMPLEControls.H"
#   include "initContinuityErrs.H"
#   include "createFields.H"
#   include "readTimeControls.H"
#   include "correctPhi.H"
#   include "CourantNo.H"
#   include "setInitialDeltaT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    scalar executionTime = runTime.elapsedCpuTime();
    scalar advectionTime = 0;

    while (runTime.run())
    {
#       include "readPIMPLEControls.H"
#       include "readTimeControls.H"
#       include "CourantNo.H"
#       include "alphaCourantNo.H"
#       include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Pressure-velocity corrector
        int oCorr = 0;
        do
        {

        scalar advectionStartTime = runTime.elapsedCpuTime();

        twoPhaseProperties.correct();

//#           include "alphaEqnSubCycle.H"

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
            alpha1 = alpha1*pos0(alpha1-clipAlphaTol)*neg0(alpha1-(1.0-clipAlphaTol)) + pos0(alpha1-(1.0-clipAlphaTol));
        }
        if ( snapAlpha )
        {
            alpha1 = min(1.0,max(0.0,alpha1));
        }

        rho == alpha1*rho1 + (scalar(1) - alpha1)*rho2;
        rhoPhi = advector.getRhoPhi(rho1, rho2);

        advectionTime += (advectionStartTime - runTime.elapsedCpuTime());


#           include "UEqn.H"

            // --- PISO loop
            for (int corr = 0; corr < nCorr; corr++)
            {
#               include "pEqn.H"
            }

#           include "continuityErrs.H"

            p = pd + rho*gh;

            if (pd.needReference())
            {
                p += dimensionedScalar
                (
                    "p",
                    p.dimensions(),
                    pRefValue - getRefCellValue(p, pdRefCell)
                );
            }

            turbulence->correct();
        } while (++oCorr < nOuterCorr);

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
