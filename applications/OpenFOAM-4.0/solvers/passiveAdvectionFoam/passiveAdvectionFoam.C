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
    passiveAdvectionFoam

Description
	Based on scalarTransportFoam
    Solves a transport equation for a passive scalar
	Modified to 
	1. run with dimensionless alpha1 field 
	2. run with adaptive timestep
	3. no diffusion term
    Used to 
	
Author
	Johan Roenby, DHI, all rights reserved
	
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
//#include "fvIOoptionList.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
	
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createTimeControls.H"

    #include "createFvOptions.H"

    simpleControl simple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "CourantNo.H"
    #include "alphaCourantNo.H"
    #include "setDeltaT.H"

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"
        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            solve
            (
                fvm::ddt(alpha1)
              + fvm::div(phi, alpha1)
//              - fvm::laplacian(DT, alpha1)
             ==
                fvOptions(alpha1)
            );
        }

		//Write total VOF and discrepancy from original VOF to log
		label lMin = -1, lMax = -1;
		scalar aMax = -GREAT, aMin = GREAT;
		forAll(alpha1,ci)
		{
			if ( alpha1[ci] > aMax)
            {
				aMax = alpha1[ci];
				lMax = ci;
			}
			else if ( alpha1[ci] < aMin )
			{
				aMin = alpha1[ci];
				lMin = ci;				
			}
		}
		
		const scalar V = sum(mesh.V()*alpha1).value();
		scalar t = runTime.time().value();
        Info << "t = " << t << ",\t sum(alpha*V) = " << V
             << ",\t dev = " << 100*(1.0-V/V0) << "%" 
             << ",\t 1-max(alpha1) = " << 1-aMax << " at cell " << lMax
             << ",\t min(alpha1) = " << aMin << " at cell " << lMin << endl;
		
        runTime.write();
		
		Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
		<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
		<< nl << endl;

    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
