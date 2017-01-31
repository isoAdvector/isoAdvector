/*---------------------------------------------------------------------------*\
|             isoAdvector | Copyright (C) 2016 Johan Roenby, DHI              |
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

Class
    Foam::isoCutFace

Description
    Class for cutting a face, faceI, of an fvMesh, mesh_, at its intersection
    with an isosurface defined by the mesh point values f_ and the isovalue,
    isoValue_.

    @article{Roenby160405,
        author={Roenby, Johan and Bredmose, Henrik and Jasak, Hrvoje},
        title={A computational method for sharp interface advection},
        volume={3},
        year={2016},
        doi={10.1098/rsos.160405},
        publisher = {The Royal Society},
        journal = {Royal Society Open Science}
    }

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

        //- Boolean telling if subface centre and area have been calculated
        bool subFaceCentreAndAreaIsCalculated_;


    // Private Member Functions

        void calcSubFaceCentreAndArea();

        //- Calculate cut points along edges of face with values f[pLabels]
        void calcSubFace
        (
            const pointField& points,
            const scalarField& f,
            const labelList& pLabels
        );

        void subFacePoints(const pointField& points, const labelList& pLabels);

        void surfacePoints(const pointField& points, const labelList& pLabels);


public:

    // Constructors

        //- Construct from fvMesh and a scalarField
        //  Length of scalarField should equal number of mesh points
        isoCutFace(const fvMesh& mesh, scalarField& f);


    // Member functions

        //- Calculate cut points along edges of faceI
        label calcSubFace(const label faceI, const scalar isoValue);

        //- Calculate cut points along edges of face with values f.
        //  Does not use the fvMesh or scalarField provided in the constructor.
        label calcSubFace
        (
            const pointField& points,
            const scalarField& f,
            const scalar isoValue
        );

        point subFaceCentre();

        vector subFaceArea();

        DynamicList<point> subFacePoints();

        DynamicList<point> surfacePoints();

        void clearStorage();
        
        //- Calculate time integrated area for a face
        scalar timeIntegratedArea
        (
            const pointField& fPts,
            const scalarField& pTimes,
            const scalar dt,
            const scalar magSf,
            const scalar Un0
        );
        
        DynamicList<point> cutPoints
        (
            const pointField& pts,
            const scalarField& f,
            const scalar f0
        );

        
        void quadAreaCoeffs
        (
            const DynamicList<point>& pf0,
            const DynamicList<point>& pf1,
            scalar& alpha,
            scalar& beta
        ) const;
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
