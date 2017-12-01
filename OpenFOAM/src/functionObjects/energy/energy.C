/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "energy.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(energy, 0);
    addToRunTimeSelectionTable(functionObject, energy, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::energy::energy
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    z0_
    (
        dimensionedScalar
        (
            "z0", dimLength, dict.lookupOrDefault<scalar>("z0", 0.0)
        )
    ),
    Etot0_(0.0),
    Evisc_(0.0),
    Ekin_
    (
        IOobject
        (
            "Ekin",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0.0*mesh_.lookupObject<volScalarField>("rho")
            *magSqr(mesh_.lookupObject<volVectorField>("U"))
    ),
    Epot_
    (
        IOobject
        (
            "Epot",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        -0.0*mesh_.lookupObject<volScalarField>("rho")
            *(mesh_.C() & mesh_.lookupObject<uniformDimensionedVectorField>("g"))
    )
{
    read(dict);
    write();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::energy::~energy()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::energy::read(const dictionary& dict)
{
//    dict.readIfPresent("wordData", wordData_);
//    dict.lookup("scalarData") >> scalarData_;
//    dict.lookup("labelData") >> labelData_;

    return true;
}


bool Foam::functionObjects::energy::execute()
{
    return true;
}


bool Foam::functionObjects::energy::end()
{
    return true;
}


bool Foam::functionObjects::energy::write()
{
    const volScalarField& rho =
        mesh_.lookupObject<volScalarField>("rho");

    const volVectorField& U =
        mesh_.lookupObject<volVectorField>("U");

    Ekin_ = 0.5*rho*magSqr(U);

    const uniformDimensionedVectorField& g =
         mesh_.lookupObject<uniformDimensionedVectorField>("g");

    vector zhat = -(g/mag(g)).value();
    Epot_ = -rho*(g & (mesh_.C() - z0_*zhat));

    scalar Ekin = gSum(fvc::volumeIntegrate(Ekin_));
    scalar Epot = gSum(fvc::volumeIntegrate(Epot_));
    scalar Etot = Ekin + Epot;

    const volScalarField& nu = mesh_.lookupObject<const volScalarField>("nu");
    scalar dEviscdt = gSum(fvc::volumeIntegrate(rho*nu*magSqr(fvc::grad(U))));
    scalar dEdt = (Etot - Etot0_)/mesh_.time().deltaTValue();
    Etot0_ = Etot;

    Info << "E = " << Etot << ", Ekin = " << Ekin << ", Epot = " << Epot
        << ", dEdt = " << dEdt << ", dEviscdt = " << dEviscdt << endl;


    return true;
}


// ************************************************************************* //
