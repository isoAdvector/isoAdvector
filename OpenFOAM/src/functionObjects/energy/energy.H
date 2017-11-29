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

Class
    Foam::functionObjects::energy

Group

Description
    This function object calculates the kinetic and potential energy for an
    interFoam type simulation. Energy density fields are written to time
    directories and total values are written to log.

    Example of function object specification:
    \verbatim
    energy1
    {
        type           energy;
        libs ("libenergyFunctionObject.so");
        ...
        z0             ;
    }
    \endverbatim

Usage
    \table
        Property     | Description                   | Required | Default
        type         | type name: energy             | yes      |
        z0           | vertical base for pot. energy | no       | 0
    \endtable

SourceFiles
    energy.C

Author
    Johan Roenby, Stromning, 2017 

\*---------------------------------------------------------------------------*/

#ifndef energy_H
#define energy_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class energy Declaration
\*---------------------------------------------------------------------------*/

class energy
:
    public fvMeshFunctionObject
{
    // Private data

        //- scalar
        dimensionedScalar z0_;

        //- scalar
        scalar Etot0_;

        //- Viscous energy loss
        scalar Evisc_;

        //- Kinetic energy density
        volScalarField Ekin_;

        //- Potential energy density
        volScalarField Epot_;

        //- Total energy density
//        volScalarField Etot_;

    // Private Member Functions


        //- Disallow default bitwise copy construct
        energy(const energy&);

        //- Disallow default bitwise assignment
        void operator=(const energy&);


public:

    //- Runtime type information
    TypeName("energy");


    // Constructors

        //- Construct from Time and dictionary
        energy
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~energy();


    // Member Functions

        //- Read the energy data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Execute at the final time-loop, currently does nothing
        virtual bool end();

        //- Write the energy
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
