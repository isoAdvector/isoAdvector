/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                isoAdvector | Copyright (C) 2016 Johan Roenby, DHI
-------------------------------------------------------------------------------

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
    along with IsoAdvector. If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "isoAdvection.H"

// ************************************************************************* //

template<typename Type>
Type Foam::isoAdvection::faceValue
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& f,
    const label facei
) const
{
    if (mesh_.isInternalFace(facei))
    {
        return f.primitiveField()[facei];
    }
    else
    {
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

        // Boundary face. Find out which face of which patch
        const label patchi = pbm.whichPatch(facei);

        if (patchi < 0 || patchi >= pbm.size())
        {
            FatalErrorInFunction
                << "Cannot find patch for face " << facei
                << abort(FatalError);
        }

        // Handle empty patches
        const polyPatch& pp = pbm[patchi];
        if (isA<emptyPolyPatch>(pp) || pp.empty())
        {
            return pTraits<Type>::zero;
        }

        const label patchFacei = pp.whichFace(facei);
        return f.boundaryField()[patchi][patchFacei];
    }
}


template<typename Type>
void Foam::isoAdvection::faceValue
(
    GeometricField<Type, fvsPatchField, surfaceMesh>& f,
    const label facei,
    const Type& value
) const
{
    if (mesh_.isInternalFace(facei))
    {
        f.primitiveFieldRef()[facei] = value;
    }
    else
    {
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

        // Boundary face. Find out which face of which patch
        const label patchi = pbm.whichPatch(facei);

        if (patchi < 0 || patchi >= pbm.size())
        {
            FatalErrorInFunction
                << "Cannot find patch for face " << facei
                << abort(FatalError);
        }

        // Handle empty patches
        const polyPatch& pp = pbm[patchi];
        if (isA<emptyPolyPatch>(pp) || pp.empty())
        {
            return;
        }

        const label patchFacei = pp.whichFace(facei);

        f.boundaryFieldRef()[patchi][patchFacei] = value;
    }
}


// ************************************************************************* //
