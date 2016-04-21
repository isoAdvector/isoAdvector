/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of the IsoAdvector source code library, which is an
    unofficial extension to OpenFOAM.

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
    Foam::isoCutFace

Description
    Class for cutting a face, faceI, of an fvMesh, mesh_, at its intersection 
    with an isosurface defined by the mesh point values f_ and the isovalue,
    isoValue_. 

Author
    Johan Roenby, DHI, all rights reserved.

SourceFiles
    isoCutFace.C

\*---------------------------------------------------------------------------*/

#ifndef isoCutFace_H
#define isoCutFace_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class isoCutFaces Declaration
\*---------------------------------------------------------------------------*/

class isoCutFace
{

    // Private data

        //- Mesh whose cells and faces to cut at their intersection with an isoface
        const fvMesh& mesh_;

        //- Face
        label faceI_;

        //- Isofunction values at mesh points. f_size() = mesh_.nPoints().
        scalarField& f_;

        //- Isovalue used to cut face or cell
        scalar isoValue_;

        //- Point along first cut edge where isosurface cuts edge
        scalar firstEdgeCut_;

        //- Point along last cut edge where isosurface cuts edge
        scalar lastEdgeCut_;

        //- Index in mesh_.faces()[faceI_] of first fully submerged (f > f0) face point.
        label firstFullySubmergedPoint_;

        //- Index in mesh_.faces()[faceI_] of last fully submerged (f > f0) face point.
        label nFullySubmergedPoints_;

        //- Storage for centre of subface
        point subFaceCentre_;

        //- Storage for area vector of subface
        vector subFaceArea_;

        //- Storage for "area of fluid" for subface
        scalar AreaOfFluid_;

        //- Storage for subFacePoints
        DynamicList<point> subFacePoints_;

        //- Storage for subFacePoints
        DynamicList<point> surfacePoints_;

        //- A face status label taking one of the values:
        //
        //  - -1: face is fully below the isosurface
        //  -  0: face is cut, i.e. has values larger and smaller than isoValue_
        //  - +1: face is fully above the isosurface 
        label faceStatus_;

    // Private Member Functions

    void calcSubFaceCentreAndArea();

public:

    // Constructors

        //- Construct from fvMesh and a scalarField
        //  Length of scalarField should equal number of mesh points
        isoCutFace
        (
            const fvMesh&,
            scalarField&
        );


        //Member functions

        label calcSubFace
        (
            const label faceI,
            const scalar isoValue
        );

        point subFaceCentre();

        vector subFaceArea();

        DynamicList<point> subFacePoints();

        DynamicList<point> surfacePoints();

        scalar areaOfFluid();

        void clearStorage();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
