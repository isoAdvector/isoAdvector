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
    isoAdvector

Description
    Advects a volume of fluid across an FVM mesh by fluxing fluid through its
    faces. Fluid transport across faces during a time step is estimated from
    the cell cutting of isosurfaces of the VOF field.

Author
    Johan Roenby, DHI, all rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "isoAdvection.H"
#include "isoCutter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    #include "createFields.H"
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
	#include "readIsoAdvectorControls.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	scalar V0 = sum(mesh.V()*alpha1).value();	
	
    Info<< "\nStarting time loop\n" << endl;
    isoAdvection advector(alpha1,phi,U,boundAlpha,vof2IsoTol,surfCellTol,writeToLog);

    while (runTime.run())
    {
        #include "readTimeControls.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        scalar dt = runTime.deltaT().value();
        advector.advect(dt);

        Info << "sum(alpha*V) = " << sum(mesh.V()*alpha1).value()
			 << ",\t dev = " << 100*(V0-sum(mesh.V()*alpha1).value())/V0 << "%" 
             << ",\t max(alpha1)-1 = " << max(alpha1).value()-1
             << ",\t min(alpha1) = " << min(alpha1).value() << endl;
		
		if ( clipAlphaTol > 0.0 )
		{
			alpha1 = alpha1*pos(alpha1-clipAlphaTol)*neg(alpha1-(1.0-clipAlphaTol)) + pos(alpha1-(1.0-clipAlphaTol));
		}
		if ( snapAlpha )
		{
			alpha1 = min(1.0,max(0.0,alpha1));
		}
		
        runTime.write();
		
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
