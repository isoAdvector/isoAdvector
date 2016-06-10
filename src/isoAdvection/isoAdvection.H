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
    Foam::isoAdvection

Description
    Calculates the new VOF (alpha) field after time step dt given the initial
    VOF field and a velocityField U and face fluxes phi. The fluid transport
    calculation is based on an idea of using isosurfaces to estimate the
    internal distribution of fluid in cells and advecting such isosurfaces
    across the mesh faces with the velocity field interpolated to the
    isosurfaces.

Author
    Johan Roenby, DHI, all rights reserved.

SourceFiles
    isoAdvection.C

\*---------------------------------------------------------------------------*/

#ifndef isoAdvection_H
#define isoAdvection_H

#include "isoCutCell.H"
#include "isoCutFace.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "OFstream.H"
#include "interpolationCellPoint.H"
#include "cellSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class isoAdvection
{

    // Private data

    //- fvMesh
    const fvMesh& mesh_;

    //- VOF field
    volScalarField& alpha1_;

    //- VOF field
    scalarField& alpha1In_;

    //- volPointInterpolation from cell centres to points
    volPointInterpolation vpi_;

    //- VOF field interpolated to mesh points
    scalarField ap_;

    //- Face fluxes
    const surfaceScalarField& phi_;

    //- Velocity field
    const volVectorField& U_;

    //- How many iterations of the alpha bounding function
    const label nAlphaBounds_;

    //- Tolerance for search of isoFace giving specified VOF value
    const scalar vof2IsoTol_;

    //- Cells with surfCellTol_ < alpha1 < 1 - surfCellTol_ are treated as 
    //  surface cells
    const scalar surfCellTol_;

    isoCutCell isoCutCell_;
    
    isoCutFace isoCutFace_;

    //- List of processor patch labels
    labelList procPatchLabels_;

    //- For each patch if it is a processor patch this is a list of the face labels
    //  on this patch that are downwind to a surface cell.
    //  For non-processor patches the list will be empty.
    List< DynamicList<label> > surfaceCellFacesOnProcPatches_;
    
    DynamicList<label> surfCells_;

    //- Boolean saying if a cell has been touched by the bounding step  
    boolList cellIsBounded_;

    //- True for all surface cells and their neighbours
    boolList checkBounding_;
    
    //- Storage for boundary faces downwind to a surface cell
    DynamicList<label> bsFaces_;
    
    //- Storage for 
    DynamicList<vector> bsx0_;
    
    //- Storage for 
    DynamicList<vector> bsn0_;
    
    //- Storage for 
    DynamicList<scalar> bsUn0_;

    //- Storage for 
    DynamicList<scalar> bsf0_;

    //- Smallest face area
    scalar minMagSf_;
    
    // Private Member Functions

    //- For each face calculate volumetric face transport during dt
    void timeIntegratedFlux
    (
        const scalar dt,
        surfaceScalarField& dVf
    );

    //- For face fLabel calculate volumetric face transport during dt given the isoFace data provided as input
    scalar timeIntegratedFlux
    (
        const label fLabel,
        const vector& x0,
        const vector& n0,
        const scalar Un0,
        const scalar f0,
        const scalar dt,
        const scalar phi,
        const scalar magSf
    );

    scalar timeIntegratedArea
    (
        const pointField& fPts,
        const scalarField& pTimes,
        const scalar dt,
        const scalar magSf,
        const scalar Un0
    );


    //- For cell ci return labels of faces fluxing out of this cell (based on sign of phi)
    void getDownwindFaces
    (
        const label ci,
        DynamicList<label>& downwindFaces
    );


    //Return true if faceI is downwind to cellI (baseds on phi_)
    bool isADownwindFace
    (
        const label faceI,
        const label cellI
    );


    //- For face with vertices p calculate its area and integrated area between isocutting lines with isovalues f0 and f1 given vertex values f
    void quadAreaCoeffs
    (
        const DynamicList<point>& pf0,
        const DynamicList<point>& pf1,
        scalar& quadArea,
        scalar& intQuadArea
    );

    //- Find min and max of subset of scalarField f defined by labelList label
    void subSetExtrema
    (
        const scalarField& f,
        const labelList& labels,
        scalar& fMin,
        scalar& fMax
    );


    void limitFluxes
    (
        surfaceScalarField& dVfa,
        const scalar dt
    );


    void boundFromAbove
    (
        const scalarField& alpha1,
        const scalar dt,
        surfaceScalarField& dVfcorrected,
        DynamicList<label>& correctedFaces
    );


    //- Given the face volume transport dVf calculates the total volume leaving cell cLabel
    scalar netFlux
    (
        const surfaceScalarField& dVf,
        const label cLabel
    );


    //- Determine if a cell is a surface cell
    bool isASurfaceCell
    (
        const label cellI
    );

    
    //- Return face value of a surfaceScalarField
    scalar faceValue
    (
        const  surfaceScalarField& f,
        const label fLabel
    );


    //- Set face value of a surfaceScalarField
    void faceValue
    (
        surfaceScalarField& f,
        const label fLabel,
        const scalar value
    );


    //- Return face value of a surfaceScalarField
    vector faceValue
    (
        const  surfaceVectorField& f,
        const label fLabel
    );


    //- Set face value of a surfaceScalarField
    void faceValue
    (
        surfaceVectorField& f,
        const label fLabel,
        const vector value
    );


    //- Syncronize dVf across processor boundaries using upwind value
    void syncProcPatches
    (
        surfaceScalarField& dVf,
        const surfaceScalarField& phi
    );


    void checkIfOnProcPatch
    (
        const label faceI
    );


    public:

    //- Constructors

    isoAdvection
    (
        volScalarField& alpha1,
        const surfaceScalarField& phi,
        const volVectorField& U,
        const dictionary& dict
    );


    void getTransportedVolume
    (
        const scalar dt,
        surfaceScalarField& dVf
    );
    
    void getSurfaceCells
    (
        cellSet& surfCells
    );

    void getBoundedCells
    (
        cellSet& boundCells
    );

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
