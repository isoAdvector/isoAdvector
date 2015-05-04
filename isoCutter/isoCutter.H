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
    Foam::isoCutter

Description
    Cuts all fvMesh edges that have one vertex value above and one below the
    isovalue. Cut point is determined by linear interpolation along edge.
    For each cut cell the class forms an internal face coinciding with the
    isosurface and dividing the cell into two adjoint subcells, one above and
    one below the isosurface. Furthermore the function can return the fraction
    of the cut cell volume which is below the isosurface in a volScalarField
    (below is defined as the subcell where the function value is higher than
    the isovalue). Also the area fraction below the isosurface of each face
    can be returned in a surfaceScalarField.

Author
    Johan Roenby, DHI

SourceFiles
    isoCutter.C

\*---------------------------------------------------------------------------*/

#ifndef isoCutter_H
#define isoCutter_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class isoCutter Declaration
\*---------------------------------------------------------------------------*/

class isoCutter
{

    // Private data
	
		//- Mesh whose cells and faces to cut at their intersection with an isoface
		const fvMesh& mesh_;
		
		//- Function to take isosurface from. f_.size() should be equal to mesh_.nPoints().
		const scalarField& f_;
		
		//- Isovalue
		const scalar& f0_;

    // Private Member Functions
	
	//- 
    void vofCutCell
    (
        const label& cellI,
		const scalar& vof,
		const scalar& tol,
		const label& maxIter,
		scalar& vofout,
		vector& subCellCtr
    );
	
	//- Cut a cell and represent cut by the cut face labels, and for each of these, 
	//  the cut edges and the point along the edge at which it is cut.
	//  If cell is not cut these lists will be empty.
	//	If cell is strictly below isosurface (i.e. all vertex function values are strictly above f0_) fullySubmerged will be true.
    void isoCutCell
    (
        const label&,
		const scalar&,
        DynamicList<label>&,
        DynamicList<label>&,
        DynamicList<scalar>&,
        bool&
    );

    //- For a cell search through its faces and find first cut edge
    void findACutFaceEdgePair
    (
        const label& cellI,
		const scalar&,
        label& cutFace,
        label& cutEdge,
        scalar& cutPoint,
        bool& fullySubmerged
    );

	//- Starting from a specified cell face-edge pair find next cut face-edge pair
    void findNextCut
    (
        const label& cellI,
		const scalar&,
        label& cellFace,
        label& edge,
        scalar& cutPoint
    );

	//- For a given cell face-edge pair find the other cell face-edge of the same cell edge
    void otherEdgeFace
    (
        const label& cellI,
        label& cellFace,
        label& edge
    );

	//- Construct "isoFace" point list representing the union of a cell with the isosurface
    void getIsoFace
    (
        const DynamicList<label>& cutFaces,
        const DynamicList<label>& cutEdges,
        const DynamicList<scalar>& cutPoints,
        pointField& isoPoints
    );

	//- Cuts a face at intersection with isosurface and returns lower subface
    bool getSubFace
    (
        const label& faceLabel,
		const scalar& f0,
        pointField& partSubFacePts
    );

	//- Return labels of all cell faces with all points strictly below isosurface
    void fullySubmergedFaces
    (
        const label& cellI,
		const scalar& f0,
        DynamicList<label>& fullSubFaces
    );

	
	//- Write face points to .ply file e.g. for visualisation in Paraview.
	//	No attempt is made to sort out non-unique points 
    void writeFacesToPlyFile
    (
        const List<pointField>& pfl,
        const word& fileName,
        const word& fileDir
    );

	
	//- Calculate face centres and areas from pointFields assumed to represent oriented faces
    void makeFaceCentresAndAreas
    (
        const List<pointField>& pfList,
        vectorField& fCtrs,
        vectorField& fAreas
    );

	//- Calculate face centre and area from pointField assumed to represent an oriented face
	//  Taken from primitiveMesh
	void makeFaceCentreAndArea
    (
        const pointField& p,
        vector& fCtr,
        vector& fArea
    );

	//- Read face centres and face area vectors from mesh
	void getFaceCentresAndAreas
    (
        const labelList& faceLabels,
        vectorField& fCtrs,
        vectorField& fAreas
    );

	
	//- Calculate cell centre and volume from cell face centres and area vectors.
	//  Taken from primitiveMesh
    void makeCellCentreAndVol
    (
        const vectorField& fCtrs,
        const vectorField& fAreas,
        vector& cellCtr,
        scalar& cellVol
    );

	//- Print a pointField to terminal/log file
	//  Used to visually inspect cut cell faces in matlab
    void printPoints
    (
        const pointField& pf
    );

	
	//- Print a list of pointFields to terminal/log file
	//  Used to visually inspect cut cell faces in matlab
    void printPoints
    (
        const List<pointField>& pfl
    );


	//- For a list of face labels return their pointField representations
	void getFacePoints
	(
		const labelList& faceLabels,
		DynamicList<pointField>&
	);

	
	//- For a given face label return its pointField representation
	void getFacePoints
	(
		const label,
		pointField&
	);
	

public:

    // Constructors

        //- Construct from fvMesh, scalarField and an isovalue
        //  Length of scalarField should equal number of mesh points

        isoCutter
        (
            const fvMesh&,
            const scalarField&,
            const scalar&
        );

		//Member functions

        //- Cell volume fractions beneath isosurface
        void subCellFractions
		(
			const scalar&,
			volScalarField&
		);
		
		void subCellFraction
		(
			const label&,
			const scalar&,
			scalar&
		);

		void subCellFraction
		(
			const label&,
			const scalar&,
			scalar&,
			vector&
		);

		void vofCutCells
		(
			const volScalarField& alpha1,
			const scalar& tol,
			const label& maxIter,
			volScalarField& f0
		);

        //- Faces area fractions beneath isosurface
        void subFaceFractions
		(
			const scalar& f0,
			surfaceScalarField&
		);

		void updateAlpha
		(
			const volScalarField& alpha,
			const surfaceScalarField& phi,
			const volVectorField& U,
			const scalar& dt,
			scalarField& dV
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
		
		void isoFaceCentreAndArea
		(
			const label& cLabel,
			const scalar& f0,
			vector& faceCentre,
			vector& faceArea
		);
		
		//- Writes cut cell faces to ply files and to screen
        void write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
