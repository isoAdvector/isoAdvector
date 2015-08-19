/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is not part of OpenFOAM.

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

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "OFstream.H"

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

    //- VOF field interpolated to mesh points
    scalarField ap_;

    //- Face fluxes
    const surfaceScalarField& phi_;

    //- Velocity field
    const volVectorField& U_;

    //- For each cell weather it qualifies as a surface cell
    boolList isSurfaceCell_;

    //- Switch to choose if alpha should be bounded
    const bool boundAlpha_;

    //- Tolerance for search of isoFace giving specified VOF value
    const scalar vof2IsoTol_;

    //- Cells with surfCellTol_ < alpha1 < surfCellTol_ are treated as surface cells
    const scalar surfCellTol_;

    //- Control to switch detailed writing to log file on and off
    const bool writeToLog_;

    // Private Member Functions

    void timeIntegratedFlux
    (
        const scalar& dt,
        scalarField& dVf
    );

    scalar timeIntegratedFlux
    (
        const label& fLabel,
        const vector& x0,
        const vector& n0,
        const scalar& Un0,
        const scalar& f0,
        const scalar& dt
    );

    void getOutFluxFaces
    (
        const label& ci,
        DynamicList<label>& outFluxingFaces
    );

    label otherCell
    (
        const label& fLabel,
        const label& cLabel
    );


    void quadAreaCoeffs
    (
        const label& fLabel,
        const scalarField& f,
        const scalar& f0,
        const scalar& f1,
        scalar& A,
        scalar& B
    );

    void subSetExtrema
    (
        const scalarField& f,
        const labelList& labels,
        scalar& fMin,
        scalar& fMax
    );

    void boundAlpha
    (
        scalarField& dVf,
        const scalar& dt
    );

    scalar netFlux
    (
        const scalarField& dVf,
        const label& cLabel
    );

    void findSurfaceCells
    (
        DynamicList<label>& surfaceCells
    );

    void calcIsoFace
    (
        const label& ci,
        vector& x0,
        vector& n0,
        scalar& f0,
        scalar& Un0
    );


    public:

    isoAdvection
    (
        volScalarField& alpha1,
        const surfaceScalarField& phi,
        const volVectorField& U,
        const bool& boundAlpha,
        const scalar& vof2IsoTol,
        const scalar& surfCellTol,
        const bool& writeToLog
    );

    void advect
    (
        const scalar& dt
    );

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
