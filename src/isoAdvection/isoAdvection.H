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

#include "isoCutter.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "OFstream.H"
#include "interpolationCellPoint.H"

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
	
	//- volPointInterpolation from cell centres to points
	volPointInterpolation vpi_;
	
    //- VOF field interpolated to mesh points
    scalarField ap_;

    //- Face fluxes
    const surfaceScalarField& phi_;

    //- Velocity field
    const volVectorField& U_;

    //- For each cell weather it qualifies as a surface cell
    boolList isSurfaceCell_;

    //- How many iterations of the alpha bounding function
    const label nAlphaBounds_;

    //- Tolerance for search of isoFace giving specified VOF value
    const scalar vof2IsoTol_;

    //- Cells with surfCellTol_ < alpha1 < surfCellTol_ are treated as surface cells
    const scalar surfCellTol_;

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
        const scalar dt
    );

	scalar timeIntegratedArea
	(
		const label fLabel,
		const pointField& fPts,
		const scalarField& pTimes,
		const scalar dt,
		const scalar Un0
	);
	
    //- For cell ci return labels of faces fluxing out of this cell (based on sign of phi)
    void getDownwindFaces
    (
        const label ci,
        DynamicList<label>& downwindFaces
    );


    //- For face with vertices p calculate its area and integrated area between isocutting lines with isovalues f0 and f1 given vertex values f
    void getQuadArea
    (
		const label fLabel,
        const pointField& p,
        const scalarField& f,
        const scalar f0,
        const scalar f1,
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
		const volScalarField& alpha1,
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

    //- Find all surface cells
    void findSurfaceCells
    (
        DynamicList<label>& surfaceCells
    );

    //- For cell ci calculate the isoFace centre x0, normal n0, and isovalue f0, of the isosurface cutting it at to the right volume fractions. 
    void calcIsoFace
    (
        const label ci,
        vector& x0,
        vector& n0,
        scalar& f0
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

    public:

    //- Constructor
    isoAdvection
    (
        volScalarField& alpha1,
        const surfaceScalarField& phi,
        const volVectorField& U,
        const dictionary& dict
    );

    //- Update the alpha1 field to the next time step
    void advect
    (
        const scalar dt
    );

	void getTransportedVolume
    (
        const scalar dt,
        surfaceScalarField& dVf
    );
	
	void getxSnSdotnF
    (
        surfaceVectorField& xSnSdotnF
    );
	
	void getIsoCentreAndNormal
    (
        volVectorField& Ci,
        volVectorField& Si
    );

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
