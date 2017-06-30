/*---------------------------------------------------------------------------*\
|             isoAdvector | Copyright (C) 2016-2017 DHI                       |
-------------------------------------------------------------------------------

License
    This file is part of isoAdvector, which is an extension to OpenFOAM.

    IsoAdvector is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    IsoAdvector is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with IsoAdvector. If not, see <http://www.gnu.org/licenses/>.

Application
    calcAdvectErrors

Description
    Compares VOF field at a given time to an exact VOF solution for a plane,
    cylinder or sphere.

Author
    Johan Roenby, DHI, all rights reserved.

\*---------------------------------------------------------------------------*/

#include "timeSelector.H"
#include "fvCFD.H"
#include "isoCutCell.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    Foam::argList::noBanner();
    #include "setRootCase.H"

    // Reading user specified time(s)
    Foam::Time runTime(Foam::Time::controlDictName, args);
    instantList timeDirs = timeSelector::select0(runTime, args);

    // In case no times are specified we use latest time
    if (timeDirs.size() < 1)
    {
        instantList instList = runTime.times();
        timeDirs.append(instList.last());
    }

    // Reading mesh
    Foam::fvMesh mesh
    (
        Foam::IOobject
        (
            Foam::fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );

    // Reading name of alpha field from setAlphaFieldDict
    IOdictionary dict
    (
        IOobject
        (
            "setAlphaFieldDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    const word fieldName(dict.lookup("field"));

    // Setting time to user specified time
    runTime.setTime(timeDirs[0], 0);

    // Reading alpha1 field for user specified time
    volScalarField alpha1
    (
        IOobject
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    // Getting internal alpha1 field
    const scalarField& VOF_calc = alpha1.primitiveField();

    // Instantiating alpha1 field for comparison
    volScalarField alpha1_true(alpha1);
    alpha1_true = 0;

    if (timeDirs.size() == 1)
    {
        // We assume the velocity field to be constant in space and time,
        // propagate the shape centre defined in setAlphaFieldDict to the
        // specified time and compare.

        //Setting time to first time which is not 'constant' for reading U
        instantList instList = runTime.times();
        runTime.setTime(instList[1],0);
        const scalar firstTime(runTime.time().value());

        // Reading U from first time directory (field assumed to be constant)
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
        const Foam::vector U0 = U.primitiveField()[0];
        Info << "Assuming constant velocity U = " << U0 << endl;

        const vector centre0(dict.lookup("centre"));

        // Resetting time to user specified time
        runTime.setTime(timeDirs[0], 0);

        // Calculating new surface shape centre at useer specified time
        const vector centre(centre0 + U0*(runTime.time().value() - firstTime));

        const word surfType(dict.lookup("type"));
        Info<< "Processing type '" << surfType << "'" << endl;

        // Instantiating isoCutting point values and isoValue
        scalar f0 = 0.0;
        scalarField f(0, mesh.nPoints());
        
        if (surfType == "plane")
        {
            const vector direction(dict.lookup("direction"));
            f = -(mesh.points() - centre) & (direction/mag(direction));
            f0 = 0.0;
        }
        else if (surfType == "sphere")
        {
            const scalar radius(readScalar(dict.lookup("radius")));

            f = -mag(mesh.points() - centre);
            f0 = -radius;
        }
        else if (surfType == "cylinder")
        {
            const scalar radius(readScalar(dict.lookup("radius")));
            const vector direction(dict.lookup("direction"));

            f = -sqrt
            (
                sqr(mag(mesh.points() - centre))
              - sqr(mag((mesh.points() - centre) & direction))
            );
            f0 = -radius;
        }
        else
        {
            Info << "Invalid surface type specified" << endl;
            Info << "Must be plane, cylinder or sphere." << endl;
            Info << "Aborting..." << endl;
            return 0;
        }

        // Calculating VOF_true from f = f0 isosurface
        isoCutCell icc(mesh, f);
        icc.volumeOfFluid(alpha1_true, f0);
    }
    else if (timeDirs.size() == 2)
    {
        // We assume the user wants to compare the alpha field at the two
        // specified times

        // Setting time to second user specified time
        runTime.setTime(timeDirs[1], 1); //Should it be 1 the index??????

        // Reading alpha1 field for user specified time
        volScalarField alpha2
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        alpha1_true = alpha2;

    }
    else
    {
        // Aborting
        Info << "Error: " << timeDirs.size() << " times specified." << endl;
        Info << "You must specify at most 2 times for the comparison." << endl;
        Info << "Aborting..." << endl;
        return 0;
    }

    const scalarField& VOF_true = alpha1_true.primitiveField();

    // Calculating error measures
    const scalar V_calc(sum(VOF_calc*mesh.cellVolumes()));
    const scalar V_true(sum(VOF_true*mesh.cellVolumes()));
    const scalar E1(sum(mag(VOF_calc - VOF_true)*mesh.cellVolumes()));
    const scalar E1rel(E1/V_calc);
    const scalar dV(V_calc - V_true);
    const scalar dVrel((V_calc - V_true)/(V_calc + SMALL));
    const scalar aMin(min(VOF_calc));
    const scalar aMaxMinus1(max(VOF_calc) - 1);

    // Printing error measures
    Info << "Advection errors for time: " << runTime.time().value() << endl;

    Info<< "E1       = " << E1 << endl;
    Info<< "E1rel    = " << E1rel << endl;
    Info<< "dV       = " << dV << endl;
    Info<< "dVrel    = " << dVrel << endl;
    Info<< "aMin     = " << aMin << endl;
    Info<< "aMax - 1 = " << aMaxMinus1 << endl;

    return 0;
}


// ************************************************************************* //
