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
    Foam::isoAdvection

Description
    Calculates the new VOF (alpha1) field after time step dt given the initial
    VOF field and a velocityField U and face fluxes phi. The fluid transport
    calculation is based on an idea of using isosurfaces to estimate the
    internal distribution of fluid in cells and advecting such isosurfaces
    across the mesh faces with the velocity field interpolated to the
    isosurfaces.

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
    isoAdvection.C
    isoAdvectionTemplates.C

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
#include "className.H"


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

        //- Reference to face owners
        const labelList& owner_;

        //- Reference to face neighbours
        const labelList& neighbour_;

        //- Reference to cells
        const cellList& cells_;

        //- Reference to cell volumes
        const scalarField& cellVolumes_;
        
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


        // Point interpolation data

            //- Interpolation object from cell centres to points
            volPointInterpolation vpi_;

            //- VOF field interpolated to mesh points
            scalarField ap_;


        // Switches and tolerances. Tolerances need to go into toleranceSwitches

            //- Number of alpha bounding steps
            label nAlphaBounds_;

            //- Tolerance for search of isoFace giving specified VOF value
            scalar vof2IsoTol_;

            //- Tolerance for marking of surface cells:
            //  Those with surfCellTol_ < alpha1 < 1 - surfCellTol_
            scalar surfCellTol_;
            
            //- Print isofaces in a <case>/isoFaces/isoFaces_#N.ply files.
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

            //- Smallest face area
            scalar minMagSf_;


        // Additional data for parallel runs

            //- List of processor patch labels
            DynamicLabelList procPatchLabels_;

            //- For each patch if it is a processor patch this is a list of the
            //  face labels on this patch that are downwind to a surface cell.
            //  For non-processor patches the list will be empty.
            List<DynamicLabelList> surfaceCellFacesOnProcPatches_;


    // Private Member Functions

        // Copy control

        //- Disallow default bitwise copy construct
        isoAdvection(const isoAdvection&);

        //- Disallow default bitwise copy assignment
        void operator=(const isoAdvection&);


        // Advection functions

            //- For each face calculate volumetric face transport during dt
            void timeIntegratedFlux();

            //- Calculate volumetric face transport during dt given the isoFace
            //  data provided as input for face facei
            scalar timeIntegratedFlux
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
            void getDownwindFaces
            (
                const label ci,
                DynamicLabelList& downwindFaces
            ) const;

/*
            //- Find min and max of subset of scalarField f defined by
            //  labelList label
            void subSetExtrema
            (
                const scalarField& f,
                const labelList& labels,
                scalar& fMin,
                scalar& fMax
            ) const;
*/

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
            void faceValue
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
        //  multiple calls within a single time step. Note: not tested with
        //  subcycles
        void advect();


        // Access functions

            //- Return alpha field
            const volScalarField& alpha() const
            {
                return alpha1_;
            }

            //- Return cellSet of surface cells
            void getSurfaceCells(cellSet& surfCells) const;

            //- Return cellSet of bounded cells
            void getBoundedCells(cellSet& boundCells) const;

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
        
        //- Write face points to .vtk files e.g. for visualisation in Paraview
        void isoFacesToFile
        (
            const DynamicList< List<point> >& isoFacePts,
            const word fileName,
            const word fileDir
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
