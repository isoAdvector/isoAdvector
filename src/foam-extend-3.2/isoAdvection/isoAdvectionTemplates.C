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

\*---------------------------------------------------------------------------*/

#include "isoAdvection.H"

// ************************************************************************* //

template<typename Type>
Type Foam::isoAdvection::faceValue
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& f,
    const label fLabel
)
{
    if (mesh_.isInternalFace(fLabel))
    {
        return f.internalField()[fLabel];
    }
    else
    {
        // Boundary face. Find out which face of which patch
        const label patchI = mesh_.boundaryMesh().whichPatch(fLabel);

        // Handle empty patches
        if (mesh_.boundary()[patchI].size() == 0)
        {
            return 0;
        }

        if (patchI < 0 || patchI >= mesh_.boundaryMesh().size())
        {
            FatalErrorIn
            (
                "Type isoAdvection::faceValue\n"
                "(\n"
                "    const GeometricField<Type, fvsPatchField, surfaceMesh>&\n"
                "    const label\n"
                ")\n"
            )   << "Cannot find patch for face " << fLabel
                << abort(FatalError);
        }

        const label faceI =
            mesh_.boundaryMesh()[patchI].whichFace
            (
                fLabel
            );

        return f.boundaryField()[patchI][faceI];
    }
}


template<typename Type>
void Foam::isoAdvection::faceValue
(
    GeometricField<Type, fvsPatchField, surfaceMesh>& f,
    const label fLabel,
    const Type& value
)
{
    if (mesh_.isInternalFace(fLabel))
    {
        f.internalField()[fLabel] = value;
    }
    else
    {
        // Boundary face. Find out which face of which patch
        const label patchI = mesh_.boundaryMesh().whichPatch(fLabel);

        // Handle empty patches
        if (mesh_.boundary()[patchI].size() == 0)
        {
            return;
        }

        if (patchI < 0 || patchI >= mesh_.boundaryMesh().size())
        {
            FatalErrorIn
            (
                "void isoAdvection::faceValue\n"
                "(\n"
                "    const GeometricField<Type, fvsPatchField, surfaceMesh>&\n"
                "    const label\n"
                "    const Type\n"
                ")\n"
            )   << "Cannot find patch for face " << fLabel
                << abort(FatalError);
        }

        const label faceI =
            mesh_.boundaryMesh()[patchI].whichFace
            (
                fLabel
            );

        f.boundaryField()[patchI][faceI] = value;
    }
}


// ************************************************************************* //
