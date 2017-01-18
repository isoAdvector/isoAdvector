/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::vofCompressionLimiter

Description
    VOF interface compression scheme.

Author
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    vofCompression.C

\*---------------------------------------------------------------------------*/

#ifndef vofCompression_H
#define vofCompression_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class vofCompressionLimiter Declaration
\*---------------------------------------------------------------------------*/

class vofCompressionLimiter
{

public:

    vofCompressionLimiter(Istream&)
    {}

    scalar limiter
    (
        const scalar cdWeight,
        const scalar faceFlux,
        const scalar& phiP,
        const scalar& phiN,
        const vector&,
        const scalar&
    ) const
    {
        scalar pP = min(max(0, phiP), 1);
        scalar pN = min(max(0, phiN), 1);

        // Quadratic compression scheme
        // return 4*min(pP*(1 - pP), pN*(1 - pN));

        // Quartic compression scheme
        return 1 - max
            (
                sqr(1 - 4*pP*(1 - pP)),
                sqr(1 - 4*pN*(1 - pN))
            );
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
