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
    Foam::isoCutCell

Description
    Class for cutting a cell, cellI, of an fvMesh, mesh_, at its intersection 
    with an isosurface defined by the mesh point values f_ and the isovalue,
    isoValue_. 

Author
    Johan Roenby, DHI, all rights reserved.

SourceFiles
    isoCutCell.C

\*---------------------------------------------------------------------------*/

#ifndef isoCutCell_H
#define isoCutCell_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "OFstream.H"
#include "isoCutFace.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class isoCutCell Declaration
\*---------------------------------------------------------------------------*/

class isoCutCell
{

    // Private data

        //- Mesh whose cells and faces to cut at their intersection with an
        //  isosurface.
        const fvMesh& mesh_;

        //- Cell to cut
        label cellI_;

        //- Isofunction values at mesh points. f_size() = mesh_.nPoints().
        scalarField& f_;

        //- Isovalue used to cut cell
        scalar isoValue_;

        //- An isoCutFace object to get access to its face cutting functionality
        isoCutFace isoCutFace_;
        
        //- List of face labels of isoCutFaces
        DynamicList<label> isoCutFaces_;

        //- List of point lists each defining an isoCutFace
        DynamicList< DynamicList<point> > isoCutFacePoints_;

        //- List of face centres for isoCutFaces
        DynamicList<point> isoCutFaceCentres_;
        
        //- List of face area vectors for isoCutFaces
        DynamicList<vector> isoCutFaceAreas_;

        //- Storage for subFace edges belonging to isoFace
        DynamicList<DynamicList<point> > isoFaceEdges_;

        //- Points constituting the cell-isosurface intersection (isoface)
        DynamicList<point> isoFacePoints_;
        
        //- Face centre of the isoface
        point isoFaceCentre_;
        
        //- Face normal of the isoface by convention pointing from high to low
        //  values (i.e. opposite of the gradient vector).
        vector isoFaceArea_;

        //- Cell centre of the subcell of cellI which is "fully submerged", i.e.
        //  where the function value is higher than the isoValue_
        point subCellCentre_;
        
        //- Volume of fully submerged subcell
        scalar subCellVolume_;
        
        //- Volume of Fluid for cellI (subCellVolume_/mesh_.V()[cellI])
        scalar VOF_;

        //- List of fully submerged faces
        DynamicList<label> fullySubFaces_;

        //- A cell status label taking one of the values:
        //
        //  - -1: cell is fully below the isosurface
        //  -  0: cell is cut
        //  - +1: cell is fully above the isosurface 
        label cellStatus_;
        
        //- Boolean telling if subcell centre and volume have been calculated
        bool subCellCentreAndVolumeCalculated_;
        
        //- Boolean telling if isoface centre and area have been calculated
        bool isoFaceCentreAndAreaCalculated_;

        // Private Member Functions

    void calcSubCellCentreAndVolume();

    void calcIsoFaceCentreAndArea();
    
    void calcIsoFacePointsFromEdges();

public:

    // Constructors

        //- Construct from fvMesh and a scalarField
        //  Length of scalarField should equal number of mesh points
        isoCutCell
        (
            const fvMesh&,
            scalarField&
        );


        //Member functions

        label calcSubCell
        (
            const label cellI,
            const scalar isoValue
        );

        point subCellCentre();

        scalar subCellVolume();

        DynamicList<point> isoFacePoints();

        point isoFaceCentre();

        vector isoFaceArea();

        scalar VolumeOfFluid();

        scalar isoValue();

        void clearStorage();
        
        scalar vofCutCell
        (
            const label cellI,
            const scalar alpha1,
            const scalar tol,
            const label maxIter
        );

        label vofCutCell2
        (
            const label cellI,
            const scalar alpha1,
            const scalar tol,
            const label maxIter
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
