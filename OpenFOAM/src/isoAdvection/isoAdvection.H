/*---------------------------------------------------------------------------*\
              Original work | Copyright (C) 2016-2017 DHI
              Modified work | Copyright (C) 2016-2017 OpenCFD Ltd.
              Modified work | Copyright (C) 2017-2018 Johan Roenby
-------------------------------------------------------------------------------

License
    This file is part of isoAdvector which is an extension to OpenFOAM.

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
    VOF field and a velocity field U and face fluxes phi. The fluid transport
    calculation is based on an idea of using isosurfaces to estimate the
    internal distribution of fluid in cells and advecting such isosurfaces
    across the mesh faces with the velocity field interpolated to the
    isosurfaces.

    Reference:
        \verbatim
            Roenby, J., Bredmose, H. and Jasak, H. (2016).
            A computational method for sharp interface advection
            Royal Society Open Science, 3
            doi 10.1098/rsos.160405
        \endverbatim

    Original code supplied by Johan Roenby, DHI (2016)

SourceFiles
    isoAdvection.C
    isoAdvectionTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef isoAdvection_H
#define isoAdvection_H

#include "fvMesh.H"
#include "volFieldsFwd.H"
#include "surfaceFields.H"
#include "className.H"
#include "isoCutCell.H"
#include "isoCutFace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class isoAdvection
{
    // Private data types

        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;


    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Dictionary for isoAdvection controls
        const dictionary dict_;

        //- VOF field
        volScalarField& alpha1_;

        //- Often used reference to alpha1 internal field
        scalarField& alpha1In_;

        //- Reference to flux field
        const surfaceScalarField& phi_;

        //- Reference to velocity field
        const volVectorField& U_;

        //- Face volumetric water transport
        surfaceScalarField dVf_;

        //- Time spent performing interface advection
        scalar advectionTime_;



        // Point interpolation data

            //- VOF field interpolated to mesh points
            scalarField ap_;


        // Switches and tolerances. Tolerances need to go into toleranceSwitches

            //- Number of alpha bounding steps
            label nAlphaBounds_;

            //- Tolerance for search of isoFace giving specified VOF value
            scalar isoFaceTol_;

            //- Tolerance for marking of surface cells:
            //  Those with surfCellTol_ < alpha1 < 1 - surfCellTol_
            scalar surfCellTol_;

            //- Switch controlling whether to use isoface normals for interface
            //  orientation (default corresponding to false) to base it on
            //  a smoothed gradient of alpha calculation (giving better results
            //  on tri on tet meshes).
            bool gradAlphaBasedNormal_;

            //- Print isofaces in a <case>/isoFaces/isoFaces_#N.vtk files.
            //  Intended for debugging
            bool writeIsoFacesToFile_;

        // Cell and face cutting

            //- List of surface cells
            DynamicLabelList surfCells_;

            //- Cell cutting object
            isoCutCell isoCutCell_;

            //- Face cutting object
            isoCutFace isoCutFace_;

            //- Bool list for cells that have been touched by the bounding step
            boolList cellIsBounded_;

            //- True for all surface cells and their neighbours
            boolList checkBounding_;

            //- Storage for boundary faces downwind to a surface cell
            DynamicLabelList bsFaces_;

            //- Storage for boundary surface iso face centre
            DynamicVectorList bsx0_;

            //- Storage for boundary surface iso face normal
            DynamicVectorList bsn0_;

            //- Storage for boundary surface iso face speed
            DynamicScalarList bsUn0_;

            //- Storage for boundary surface iso value
            DynamicScalarList bsf0_;

        // Additional data for parallel runs

            //- List of processor patch labels
            DynamicLabelList procPatchLabels_;

            //- For each patch if it is a processor patch this is a list of the
            //  face labels on this patch that are downwind to a surface cell.
            //  For non-processor patches the list will be empty.
            List<DynamicLabelList> surfaceCellFacesOnProcPatches_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        isoAdvection(const isoAdvection&);

        //- Disallow default bitwise copy assignment
        void operator=(const isoAdvection&);


        // Advection functions

            //- For each face calculate volumetric face transport during dt
            void timeIntegratedFlux();

            //- Calculate volumetric face transport during dt given the isoFace
            //  data provided as input for face facei
            scalar timeIntegratedFaceFlux
            (
                const label facei,
                const vector& x0,
                const vector& n0,
                const scalar Un0,
                const scalar f0,
                const scalar dt,
                const scalar phi,
                const scalar magSf
            );

            //- For a given cell return labels of faces fluxing out of this cell
            //  (based on sign of phi)
            void setDownwindFaces
            (
                const label celli,
                DynamicLabelList& downwindFaces
            ) const;

            // Limit fluxes
            void limitFluxes();

            // Bound fluxes
            void boundFromAbove
            (
                const scalarField& alpha1,
                surfaceScalarField& dVfcorrected,
                DynamicLabelList& correctedFaces
            );

            //- Given the face volume transport dVf calculates the total volume
            //  leaving a given cell. Note: cannot use dVf member because
            //  netFlux is called also for corrected dVf
            scalar netFlux
            (
                const surfaceScalarField& dVf,
                const label celli
            ) const;

            //- Determine if a cell is a surface cell
            bool isASurfaceCell(const label celli) const
            {
                return
                (
                    surfCellTol_ < alpha1In_[celli]
                 && alpha1In_[celli] < 1 - surfCellTol_
                );
            }

            //- Clear out isoFace data
            void clearIsoFaceData()
            {
                surfCells_.clear();
                bsFaces_.clear();
                bsx0_.clear();
                bsn0_.clear();
                bsUn0_.clear();
                bsf0_.clear();
                
                if (mesh_.topoChanging())
                {
                    // Introduced resizing to cope with changing meshes
                    checkBounding_.resize(mesh_.nCells());
                    cellIsBounded_.resize(mesh_.nCells());
                }
                checkBounding_ = false;
                cellIsBounded_ = false;

            }

        // Face value functions needed for random face access where the face
        // can be either internal or boundary face

            //- Return face value for a given Geometric surface field
            template<typename Type>
            Type faceValue
            (
                const GeometricField<Type, fvsPatchField, surfaceMesh>& f,
                const label facei
            ) const;

            //- Set face value for a given Geometric surface field
            template<typename Type>
            void setFaceValue
            (
                GeometricField<Type, fvsPatchField, surfaceMesh>& f,
                const label facei,
                const Type& value
            ) const;


        // Parallel run handling functions

            //- Synchronize dVf across processor boundaries using upwind value
            void syncProcPatches
            (
                surfaceScalarField& dVf,
                const surfaceScalarField& phi
            );

            //- Check if the face is on processor patch and append it to the
            //  list of surface cell faces on processor patches
            void checkIfOnProcPatch(const label facei);


public:

    //- Runtime type information
    TypeName("isoAdvection");

    //- Constructors

        //- Construct given alpha, phi and velocity field. Note: phi should be
        // divergence free up to a sufficient tolerance
        isoAdvection
        (
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U
        );


    //- Destructor
    virtual ~isoAdvection()
    {}


    // Member functions

        //- Advect the free surface. Updates alpha field, taking into account
        //  multiple calls within a single time step.
        void advect();

        //- Apply the bounding based on user inputs
        void applyBruteForceBounding();

        // Access functions

            //- Return alpha field
            const volScalarField& alpha() const
            {
                return alpha1_;
            }

            //- Return the controls dictionary
            const dictionary& dict() const
            {
                return dict_;
            }

            //- Return cellSet of surface cells
            void writeSurfaceCells() const;

            //- Return cellSet of bounded cells
            void writeBoundedCells() const;

            //- Return mass flux
            tmp<surfaceScalarField> getRhoPhi
            (
                const dimensionedScalar rho1,
                const dimensionedScalar rho2
            ) const
            {
                return tmp<surfaceScalarField>
                (
                    new surfaceScalarField
                    (
                        "rhoPhi",
                        (rho1 - rho2)*dVf_/mesh_.time().deltaT() + rho2*phi_
                    )
                );
            }

            scalar advectionTime() const
            {
                return advectionTime_;
            }

            //- Write isoface points to .obj file
            void writeIsoFaces
            (
                const DynamicList<List<point> >& isoFacePts
            ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "isoAdvectionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
