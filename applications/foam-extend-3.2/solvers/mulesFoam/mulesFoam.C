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
    mulesFoam

Description
    This solver is just interFoam taking a number of extra input paramters.
	The sole purpose of the solver is to compare the behaviour of MULES for
	pure advection problems with the corresponding IsoAdvector results.

Author
    Johan Roenby, DHI, all rights reserved

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "MULES.H"
#include "subCycle.H"
#include "interfaceProperties.H"
#include "twoPhaseMixture.H"
#include "turbulenceModel.H"

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

    while (runTime.run())
    {
#       include "readPIMPLEControls.H"
#       include "readTimeControls.H"
#       include "CourantNo.H"
#       include "alphaCourantNo.H"
#       include "setDeltaT.H"

		scalar t = runTime.time().value();
        scalar dt = runTime.deltaT().value();
		if ( reverseTime > 0.0 && t >= reverseTime )
		{
			Info<< "Reversing flow" << endl;
			phi = -phi;
			phi0 = -phi0;
			U = -U;
			U0 = -U0;
			reverseTime = -1.0;
		}
		if ( period > 0.0 )
		{
			phi = phi0*Foam::cos(2.0*M_PI*(t + 0.5*dt)/period);
			U = U0*Foam::cos(2.0*M_PI*(t + 0.5*dt)/period);
		}

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Pressure-velocity corrector
        int oCorr = 0;
        do
        {
            twoPhaseProperties.correct();

#           include "alphaEqnSubCycle.H"

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

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
