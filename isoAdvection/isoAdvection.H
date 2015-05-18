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

	// Private Member Functions

	void timeIntegratedFlux
	(
		const volScalarField& alpha,
		const surfaceScalarField& phi,
		const volVectorField& U,
		const scalar& dt,
		scalarField& dVf
	);

	void getOutFluxFaces
	(
		const surfaceScalarField& phi,
		const label& ci,
		DynamicList<label>& outFluxingFaces
	);

	void getVertexAlphas
	(
		const scalarField& alphap,
		const fvMesh& mesh_,
		const label& fLabel,
		const scalar& f0,
		const scalar& Un0,
		scalarField& av,
		vector& xFirst,
		vector& xLast
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
		const volScalarField& alpha,
		const surfaceScalarField& phi,
		scalarField& dVf
	);

	scalar netFlux
	(
		const scalarField& dVf,
		const fvMesh& mesh_,
		const label& cLabel
	);
	
public:

	void advect
	(
		volScalarField& alpha1,
		const surfaceScalarField& phi,
		const volVectorField& U,
		const scalar& dt
	);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
