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

\*---------------------------------------------------------------------------*/

#include "isoAdvection.H"
#include "volPointInterpolation.H"
#include "interpolationCellPoint.H"
#include "fvcSurfaceIntegrate.H"
#include "upwind.H"

// * * * * * * * * * * * * * * Debugging * * * * * * * * * * * * * //

#ifndef DebugInfo
//Taken from OpenFOAM-4.0/src/OpenFOAM/db/error/messageStream.H to make code
//compile with older OF versions.
#define DebugInfo                                                              \
    if (debug) Info
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(isoAdvection, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isoAdvection::isoAdvection
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U
)
:
    // General data
    mesh_(alpha1.mesh()),
    dict_(mesh_.solutionDict().subDict("isoAdvector")),
    alpha1_(alpha1),
    alpha1In_(alpha1.oldTime().ref()), // Need to reorganise
    phi_(phi),
    U_(U),
    dVf_
    (
        IOobject
        (
            "dVf_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimVol, 0)
    ),

    // Interpolation data
    vpi_(mesh_),
    ap_(mesh_.nPoints()),

    // Tolerances and solution controls
    nAlphaBounds_(dict_.lookupOrDefault<label>("nAlphaBounds", 3)),
    vof2IsoTol_(dict_.lookupOrDefault<scalar>("vof2IsoTol", 1e-8)),
    surfCellTol_(dict_.lookupOrDefault<scalar>("surfCellTol", 1e-8)),

    // Cell cutting data
    surfCells_(label(0.2*mesh_.nCells())),
    isoCutCell_(mesh_, ap_),
    isoCutFace_(mesh_, ap_),
    cellIsBounded_(mesh_.nCells(), false),
    checkBounding_(mesh_.nCells(), false),
    bsFaces_(label(0.2*(mesh_.nFaces() - mesh_.nInternalFaces()))),
    bsx0_(bsFaces_.size()),
    bsn0_(bsFaces_.size()),
    bsUn0_(bsFaces_.size()),
    bsf0_(bsFaces_.size()),
    minMagSf_(gMin(mesh_.magSf())),

    // Parallel run data
    procPatchLabels_(mesh_.boundary().size()),
    surfaceCellFacesOnProcPatches_(0)
{
    // Prepare lists used in parallel runs
    if (Pstream::parRun())
    {
        // Force calculation of cell centres and volumes (else parallel
        // communication may crash)
        mesh_.C();

        // Get boundary mesh and resize the list for parallel comms
        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        surfaceCellFacesOnProcPatches_.resize(patches.size());

        // Append all processor patch labels to the list
        forAll(patches, patchi)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchi])
             && patches[patchi].size() > 0
            )
            {
                procPatchLabels_.append(patchi);
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::isoAdvection::timeIntegratedFlux()
{

    // Get time step
    const scalar dt = mesh_.time().deltaT().value();

    // Create interpolation object for interpolating velocity to iso face
    // centres
    interpolationCellPoint<vector> UInterp(U_);

    // For each downwind face of each surface cell we "isoadvect" to find dVf
    label nSurfaceCells = 0;

    // Clear out the data for re-use and reset list containing information
    // whether cells could possibly need bounding
    clearIsoFaceData();
    checkBounding_ = false;

    // Get necessary references
    const scalarField& phiIn = phi_.primitiveField();
    const scalarField& magSfIn = mesh_.magSf().primitiveField();
    scalarField& dVfIn = dVf_.primitiveFieldRef();

    // Get necessary mesh data
    const labelListList& CP = mesh_.cellPoints();
    const labelListList& CC = mesh_.cellCells();
    const cellList& meshCells = mesh_.cells();
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    // Loop through cells
    forAll(alpha1In_, celli)
    {
        if (isASurfaceCell(celli))
        {
            // This is a surface cell, increment the counter, append and mark
            // the cell
            nSurfaceCells++;
            surfCells_.append(celli);
            checkBounding_[celli] = true;

            DebugInfo
                << "\n------------ Cell " << celli << " with alpha1 = "
                << alpha1In_[celli] << " and 1-alpha1 = "
                << 1.0 - alpha1In_[celli] << " ------------"
                << endl;

            // Calculate isoFace centre x0, normal n0 at time t
            label maxIter = 100; // NOTE: make it a debug switch

            // Calculate cell status (-1: cell is fully below the isosurface, 0:
            // cell is cut, 1: cell is fully above the isosurface)
            label cellStatus = isoCutCell_.vofCutCell
            (
                celli,
                alpha1In_[celli],
                vof2IsoTol_,
                maxIter
            );

//            Info << "1 - f0 = " << 1 - f0 << " for cell " << celli << endl;

            // Cell is cut
            if (cellStatus == 0)
            {
                const scalar f0 = isoCutCell_.isoValue();
                const point x0 = isoCutCell_.isoFaceCentre();
                vector n0 = isoCutCell_.isoFaceArea();

                // If cell almost full or empty isoFace may be undefined.
                // Calculating normal by going a little into the cell.
                // Note: put 1e-6 into the tolerance
                if (mag(n0) < 1e-6*minMagSf_)
                {
                    WarningInFunction
                        << "mag(n0) = " << mag(n0)
                        << " < 1e-6*minMagSf_ for cell " << celli << endl;

                    // Initialise minimum and maximum values
                    scalar fMin = GREAT;
                    scalar fMax = -GREAT;

                    // Get cell points
                    const labelList& cellPts = CP[celli];

                    // Calculate min and max values of the subset
                    subSetExtrema(ap_, cellPts, fMin, fMax);

                    scalar fInside  = 0;
                    if (alpha1In_[celli] > 0.5) // Note: changed from >= to >
                    {
                        // Note: make a tolerance
                        fInside =  fMin + 1e-3*(fMax-fMin);
                    }
                    else
                    {
                        fInside =  fMax - 1e-3*(fMax-fMin);
                    }

//                    scalar fInside = f0 + sign(alpha1In_[celli]-0.5)*1e-3;

                    // Calculate sub cell and initialise the normal with face
                    // area vector
                    isoCutCell_.calcSubCell(celli,fInside);
                    n0 = isoCutCell_.isoFaceArea();
                }

                if (mag(n0) > 1e-6*minMagSf_)
                {
                    DebugInfo << "Normalising n0: " << n0 << endl;
                    n0 /= mag(n0);
                }
                else
                {
                    WarningInFunction
                        << "mag(n0) = " << mag(n0)
                        << " < 1e-6*minMagSf_ for cell " << celli
                        << " with alpha1 = " << alpha1In_[celli]
                        << ", 1-alpha1 = " << 1.0-alpha1In_[celli]
                        << " and f0 = " << f0 << endl;

                    // Normalise the vector with stabilisation
                    n0 /= (mag(n0) + SMALL);
                    Info<< "After normalisation: mag(n0) = "
                        << mag(n0) << endl;
                }

                // Get the speed of the isoface by interpolating velocity and
                // dotting it with isoface normal
                const scalar Un0 = UInterp.interpolate(x0, celli) & n0;

                DebugInfo
                    << "calcIsoFace gives initial surface: \nx0 = " << x0
                    << ", \nn0 = " << n0 << ", \nf0 = " << f0 << ", \nUn0 = "
                    << Un0 << endl;

                // Estimating time integrated water flux through each downwind
                // face
                const cell& cellFaces = meshCells[celli];
                forAll(cellFaces, fi)
                {
                    // Get current face index
                    const label facei = cellFaces[fi];

                    // Check if the face is internal face
                    if (mesh_.isInternalFace(facei))
                    {
                        bool isDownwindFace = false;
                        label otherCell = -1;

                        // Check if the cell is owner
                        if (celli == own[facei])
                        {

                            if (phiIn[facei] > 10*SMALL)
                            {
                                isDownwindFace = true;
                            }

                            // Other cell is neighbour
                            otherCell = nei[facei];
                        }
                        else //celli is the neighbour
                        {
                            if (phiIn[facei] < -10*SMALL)
                            {
                                isDownwindFace = true;
                            }

                            // Other cell is the owner
                            otherCell = own[facei];
                        }

                        // Calculate time integrated flux if this is a downwind
                        // face
                        if (isDownwindFace)
                        {
//                            Info<< "Setting value for internal face " << facei
//                                << endl;
                            dVfIn[facei] = timeIntegratedFlux
                            (
                                facei,
                                x0,
                                n0,
                                Un0,
                                f0,
                                dt,
                                phiIn[facei],
                                magSfIn[facei]
                            );
                        }

                        // We want to check bounding of neighbour cells to surface
                        // cells as well:
                        checkBounding_[otherCell] = true;

                        // Also check neighbours of neighbours.
                        // Note: consider making it a run time selectable
                        // extension level (easily done with recursion):
                        // 0 - only neighbours
                        // 1 - neighbours of neighbours
                        // 2 - ...
                        const labelList& nNeighbourCells = CC[otherCell];
                        forAll(nNeighbourCells, ni)
                        {
                            checkBounding_[nNeighbourCells[ni]] = true;
                        }
                    }
                    else
                    {
                        bsFaces_.append(facei);
                        bsx0_.append(x0);
                        bsn0_.append(n0);
                        bsUn0_.append(Un0);
                        bsf0_.append(f0);
                        // Note: we must not check if the face is on the
                        // processor patch here.
                    }
                }
            }
        }
    }

    // Get references to boundary fields
    const surfaceScalarField::Boundary& phib = phi_.boundaryField();
    const surfaceScalarField::Boundary& magSfb = mesh_.magSf().boundaryField();
    surfaceScalarField::Boundary& dVfb = dVf_.boundaryFieldRef();

    // Loop through boundary surface faces
    forAll(bsFaces_, fi)
    {
        // Get boundary face index (in the global list)
        const label facei = bsFaces_[fi];

        // Get necesary mesh data
        const fvBoundaryMesh& boundaryMesh = mesh_.boundary();
        const polyBoundaryMesh& pBoundaryMesh = mesh_.boundaryMesh();

        // Get necessary labels
        // Note: consider optimisation since whichPatch is expensive
        const label patchi = pBoundaryMesh.whichPatch(facei);
        const label start = boundaryMesh[patchi].patch().start();
        const label size = boundaryMesh[patchi].size();

        if (size > 0)
        {
            // Get patch local label
            const label patchFacei = facei - start;
            const scalar phiP = phib[patchi][patchFacei];

            if (phiP > 0) // Note: changed from phiP > 10*SMALL
            {
                const scalar magSf = magSfb[patchi][patchFacei];

                dVfb[patchi][patchFacei] = timeIntegratedFlux
                (
                    facei,
                    bsx0_[fi],
                    bsn0_[fi],
                    bsUn0_[fi],
                    bsf0_[fi],
                    dt,
                    phiP,
                    magSf
                );

                // Check if the face is on processor patch and append it to
                // the list if necessary
                checkIfOnProcPatch(facei);
            }
        }
    }

    // Print out number of surface cells
    Info<< "Number of isoAdvector surface cells = "
        << returnReduce(nSurfaceCells, sumOp<label>()) << endl;
}


Foam::scalar Foam::isoAdvection::timeIntegratedFlux
(
    const label facei,
    const vector& x0,
    const vector& n0,
    const scalar Un0,
    const scalar f0,
    const scalar dt,
    const scalar phi,
    const scalar magSf
)
{
    // Note: this function is often called within a loop. Consider passing mesh
    // faces, volumes and points as arguments instead of accessing here

    //Treating rare cases where isoface normal is not calculated properly
    if (mag(n0) < 0.5)
    {
        scalar alphaf = 0.0;
        scalar waterInUpwindCell = 0.0;

        if (phi > 0 || !mesh_.isInternalFace(facei))
        {
            const label upwindCell = mesh_.faceOwner()[facei];
            alphaf = alpha1In_[upwindCell];
            waterInUpwindCell = alphaf*mesh_.V()[upwindCell];
        }
        else
        {
            const label upwindCell = mesh_.faceNeighbour()[facei];
            alphaf = alpha1In_[upwindCell];
            waterInUpwindCell = alphaf*mesh_.V()[upwindCell];
        }

        if (debug)
        {
            WarningInFunction
                << "mag(n0) = " << mag(n0)
                << " so timeIntegratedFlux calculates dVf from upwind"
                << " cell alpha value: " << alphaf << endl;
        }

        return min(alphaf*phi*dt, waterInUpwindCell);
    }


    // Find sorted list of times where the isoFace will arrive at face points
    // given initial position x0 and velocity Un0*n0

    // Get points for this face
    const face& f = mesh_.faces()[facei];
    const pointField fPts(f.points(mesh_.points()));
    const label nPoints = fPts.size();

    scalarField pTimes(fPts.size());
    if (mag(Un0) > 1e-12) // Note: tolerances
    {
        // Here we estimate time of arrival to the face points from their normal
        // distance to the initial surface and the surface normal velocity

        pTimes = ((fPts - x0) & n0)/Un0;

        scalar dVf = 0;

        // Check if pTimes changes direction more than twice when looping face
        label nShifts = 0;
        forAll(pTimes, pi)
        {
            const label oldEdgeSign =
                sign(pTimes[(pi + 1) % nPoints] - pTimes[pi]);
            const label newEdgeSign =
                sign(pTimes[(pi + 2) % nPoints] - pTimes[(pi + 1) % nPoints]);

            if (newEdgeSign != oldEdgeSign)
            {
                nShifts++;
            }
        }
        if (nShifts == 2)
        {
            dVf = phi/magSf*isoCutFace_.timeIntegratedArea(fPts, pTimes, dt, magSf, Un0);
        }
        else if (nShifts > 2) // triangle decompose the face
        {
            pointField fPts_tri(3);
            scalarField pTimes_tri(3);
            fPts_tri[0] = mesh_.faceCentres()[facei];
            pTimes_tri[0] = ((fPts_tri[0] - x0) & n0)/Un0;
            for (label pi = 0; pi < nPoints; pi++)
            {
                fPts_tri[1] = fPts[pi];
                pTimes_tri[1] = pTimes[pi];
                fPts_tri[2] = fPts[(pi + 1) % nPoints];
                pTimes_tri[2] = pTimes[(pi + 1) % nPoints];
                const scalar magSf_tri =
                    mag
                    (
                        0.5*(fPts_tri[2] - fPts_tri[0])^(fPts_tri[1]
                      - fPts_tri[0])
                    );
                const scalar phi_tri = phi*magSf_tri/magSf;
                dVf +=
                    phi_tri
                   /magSf_tri
                   *isoCutFace_.timeIntegratedArea
                    (
                        fPts_tri,
                        pTimes_tri,
                        dt,
                        magSf_tri,
                        Un0
                    );
            }
//            WarningInFunction
//                << "Warning: nShifts = " << nShifts << " on face " << facei
//                << " owned by cell " << mesh_.faceOwner()[facei] << endl;
        }
        else
        {
            WarningInFunction
                << "Warning: nShifts = " << nShifts << " on face " << facei
                << " with pTimes = " << pTimes << " owned by cell "
                << mesh_.faceOwner()[facei] << endl;
        }

        return dVf;
    }
    else
    {
        // Un0 is almost zero and isoFace is treated as stationary
        isoCutFace_.calcSubFace(facei, f0);
        const scalar alphaf = mag(isoCutFace_.subFaceArea()/magSf);

        if (debug)
        {
            WarningInFunction
                << "Un0 is almost zero (" << Un0
                << ") - calculating dVf on face " << facei
                << " using subFaceFraction giving alphaf = " << alphaf
                << endl;
        }

        return phi*dt*alphaf;
    }
}


void Foam::isoAdvection::getDownwindFaces
(
    const label celli,
    DynamicLabelList& downwindFaces
) const
{

    // Get necessary mesh data and cell information
    const labelList& own = mesh_.faceOwner();
    const cellList& cells = mesh_.cells();
    const cell& c = cells[celli];

    // Check all faces of the cell
    forAll(c, fi)
    {
        // Get face and corresponding flux
        const label facei = c[fi];
        const scalar& phi = faceValue(phi_, facei);

        if (own[facei] == celli)
        {
            if (phi > 10*SMALL)
            {
                downwindFaces.append(facei);
            }
        }
        else if (phi < -10*SMALL)
        {
            downwindFaces.append(facei);
        }
    }

    downwindFaces.shrink();
}


void Foam::isoAdvection::subSetExtrema
(
    const scalarField& f,
    const labelList& labels,
    scalar& fMin,
    scalar& fMax
) const
{
    fMin = VGREAT;
    fMax = -VGREAT;

    forAll(labels, pi)
    {
        const scalar fp = f[labels[pi]];
        fMin = min(fMin, fp);
        fMax = max(fMax, fp);
    }
}


void Foam::isoAdvection::limitFluxes()
{

    // Get time step size
    const scalar dt = mesh_.time().deltaT().value();

//    scalarField alphaNew = alpha1In_ - fvc::surfaceIntegrate(dVf_);
    const scalar aTol = 1.0e-12;          // Note: tolerances
    const scalar maxAlphaMinus1 = 1;      // max(alphaNew - 1);
    const scalar minAlpha = -1;           // min(alphaNew);
    const label nUndershoots = 20;        // sum(neg(alphaNew + aTol));
    const label nOvershoots = 20;         // sum(pos(alphaNew - 1 - aTol));
    cellIsBounded_ = false;

    // Loop number of bounding steps
    for (label n = 0; n < nAlphaBounds_; n++)
    {
        Info<< "isoAdvection: bounding iteration " << n + 1 << endl;

        if (maxAlphaMinus1 > aTol) // Note: tolerances
        {
            DebugInfo << "Bound from above... " << endl;

//          scalarField& dVfcorrected = dVf_.primitiveFieldRef();

            surfaceScalarField dVfcorrected("dVfcorrected", dVf_);
            DynamicList<label> correctedFaces(3*nOvershoots);
            boundFromAbove(alpha1In_, dVfcorrected, correctedFaces);

            forAll(correctedFaces, fi)
            {
                label facei = correctedFaces[fi];

                // Change to treat boundaries consistently
                faceValue(dVf_, facei, faceValue(dVfcorrected, facei));
            }

            syncProcPatches(dVf_, phi_);
        }

        if (minAlpha < -aTol) // Note: tolerances
        {
            DebugInfo << "Bound from below... " << endl;

            scalarField alpha2 = 1.0 - alpha1In_;
            surfaceScalarField dVfcorrected
            (
                "dVfcorrected",
                phi_*dimensionedScalar("dt", dimTime, dt) - dVf_
            );
//          dVfcorrected -= dVf_;   // phi_ and dVf_ have same sign and dVf_ is
                                    // the portion of phi_*dt that is water.
            // If phi_ > 0 then dVf_ > 0 and mag(phi_*dt-dVf_) < mag(phi_*dt) as
            // it should.
            // If phi_ < 0 then dVf_ < 0 and mag(phi_*dt-dVf_) < mag(phi_*dt) as
            // it should.
            DynamicList<label> correctedFaces(3*nUndershoots);
            boundFromAbove(alpha2, dVfcorrected, correctedFaces);
            forAll(correctedFaces, fi)
            {
                label facei = correctedFaces[fi];

                // Change to treat boundaries consistently
                scalar phi = faceValue(phi_, facei);
                scalar dVcorr = faceValue(dVfcorrected, facei);
                faceValue(dVf_, facei, phi*dt - dVcorr);
            }

            syncProcPatches(dVf_, phi_);
        }

//        // Check if still unbounded
//        alphaNew = alpha1In_ - fvc::surfaceIntegrate(dVf);
//        maxAlphaMinus1 = max(alphaNew-1);
//        minAlpha = min(alphaNew);
//        nUndershoots = sum(neg(alphaNew+aTol));
//        nOvershoots = sum(pos(alphaNew-1-aTol));
//        Info<< "After bounding number " << n + 1 << " of time "
//            << mesh_.time().value() << ":" << endl;
//        Info<< "nOvershoots = " << nOvershoots << " with max(alphaNew-1) = "
//            << maxAlphaMinus1 << " and nUndershoots = " << nUndershoots
//            << " with min(alphaNew) = " << minAlpha << endl;
    }
}


void Foam::isoAdvection::boundFromAbove
(
    const scalarField& alpha1,
    surfaceScalarField& dVf,
    DynamicList<label>& correctedFaces
)
{

    // Get time step size
    const scalar dt = mesh_.time().deltaT().value();

    correctedFaces.clear();
    scalar aTol = 10*SMALL; // Note: tolerances
    scalar maxOvershoot = -GREAT;
    label maxOvershootCell = -1;
    
    // Get necessary mesh data
    const scalarField& meshV = mesh_.V();
    const cellList& meshCells = mesh_.cells();

    // Loop through alpha cell centred field
    forAll(alpha1, celli)
    {
        if (checkBounding_[celli])
        {
            const scalar& Vi = meshV[celli];
            scalar alpha1New = alpha1[celli] - netFlux(dVf, celli)/Vi;
            scalar alphaOvershoot = alpha1New - 1.0;
            scalar fluidToPassOn = alphaOvershoot*Vi;
            label nFacesToPassFluidThrough = 1;
            
            if (alphaOvershoot > maxOvershoot)
            {
                maxOvershoot = alphaOvershoot;
                maxOvershootCell = celli;
            }

            bool firstLoop = true;
            // First try to pass surplus fluid on to neighbour cells that are
            // not filled and to which dVf < phi*dt
            while (alphaOvershoot > aTol && nFacesToPassFluidThrough > 0)
            {
                DebugInfo
                    << "\n\nBounding cell " << celli
                    << " with alpha overshooting " << alphaOvershoot
                    << endl;

                cellIsBounded_[celli] = true;

                // Find potential neighbour cells to pass surplus phase to
                DynamicList<label> downwindFaces(meshCells[celli].size());
                getDownwindFaces(celli, downwindFaces);

                DynamicList<label> facesToPassFluidThrough(downwindFaces.size());
                DynamicList<scalar> dVfmax(downwindFaces.size());
                DynamicList<scalar> phi(downwindFaces.size());

                scalar dVftot = 0.0;
                nFacesToPassFluidThrough = 0;

                DebugInfo << "downwindFaces: " << downwindFaces << endl;

                forAll(downwindFaces, fi)
                {
                    const label facei = downwindFaces[fi];
                    const scalar phif = faceValue(phi_, facei);
                    const scalar dVff = faceValue(dVf, facei);
                    const scalar maxExtraFaceFluidTrans = mag(phif*dt - dVff);

                    // dVf has same sign as phi and so if phi>0 we have
                    // mag(phi_[facei]*dt) - mag(dVf[facei]) = phi_[facei]*dt
                    // - dVf[facei]
                    // If phi < 0 we have mag(phi_[facei]*dt) -
                    // mag(dVf[facei]) = -phi_[facei]*dt - (-dVf[facei]) > 0
                    // since mag(dVf) < phi*dt
                    DebugInfo
                        << "downwindFace " << facei
                        << " has maxExtraFaceFluidTrans = "
                        << maxExtraFaceFluidTrans << endl;

                    if (maxExtraFaceFluidTrans/Vi > aTol)
//                    if (maxExtraFaceFluidTrans/Vi > aTol &&
//                    mag(dVfIn[facei])/Vi > aTol) //Last condition may be
//                    important because without this we will flux through uncut
//                    downwind faces
                    {
                        facesToPassFluidThrough.append(facei);
                        phi.append(phif);
                        dVfmax.append(maxExtraFaceFluidTrans);
                        dVftot += mag(phif*dt);
                    }
                }

                DebugInfo
                    << "\nfacesToPassFluidThrough: "
                    << facesToPassFluidThrough << ", dVftot = "
                    << dVftot << " m3 corresponding to dalpha = "
                    << dVftot/Vi << endl;

                forAll(facesToPassFluidThrough, fi)
                {
                    const label facei = facesToPassFluidThrough[fi];
                    scalar fluidToPassThroughFace =
                        fluidToPassOn*mag(phi[fi]*dt)/dVftot;

                    nFacesToPassFluidThrough +=
                        pos(dVfmax[fi] - fluidToPassThroughFace);

                    fluidToPassThroughFace =
                        min(fluidToPassThroughFace, dVfmax[fi]);

                    scalar dVff = faceValue(dVf, facei);
                    dVff += sign(phi[fi])*fluidToPassThroughFace;
                    faceValue(dVf, facei, dVff);

                    if (firstLoop)
                    {
                        checkIfOnProcPatch(facei);
                        correctedFaces.append(facei);
                    }
                }

                firstLoop = false;
                alpha1New = alpha1[celli] - netFlux(dVf, celli)/Vi;
                alphaOvershoot = alpha1New - 1.0;
                fluidToPassOn = alphaOvershoot*Vi;

                DebugInfo
                    << "\nNew alpha for cell " << celli << ": "
                    << alpha1New << endl;
            }
        }
    }

//    Info << "maxOvershoot = " << maxOvershoot << endl;
    DebugInfo << "correctedFaces = " << correctedFaces << endl;
}


Foam::scalar Foam::isoAdvection::netFlux
(
    const surfaceScalarField& dVf,
    const label celli
) const
{
    scalar dV = 0.0;

    // Get face label
    const cell& c = mesh_.cells()[celli];

    // Get mesh data
    const labelList& own = mesh_.faceOwner();

    forAll(c, fi)
    {
        const label facei = c[fi];
        const scalar dVff = faceValue(dVf, facei);

        if (own[facei] == celli)
        {
            dV += dVff;
        }
        else
        {
            dV -= dVff;
        }
    }

    return dV;
}


void Foam::isoAdvection::syncProcPatches
(
    surfaceScalarField& dVf,
    const surfaceScalarField& phi
)
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    if (Pstream::parRun())
    {
        PstreamBuffers pBufs(Pstream::nonBlocking);

        // Send
        forAll(procPatchLabels_, patchLabeli)
        {
            const label patchi = procPatchLabels_[patchLabeli];

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchi]);

            UOPstream toNbr(procPatch.neighbProcNo(), pBufs);
            const scalarField& pFlux = dVf.boundaryField()[patchi];

            const List<label>& surfCellFacesOnProcPatch =
                surfaceCellFacesOnProcPatches_[patchi];

            List<scalar> dVfPatch(surfCellFacesOnProcPatch.size());
            forAll(dVfPatch, facei)
            {
                dVfPatch[facei] = pFlux[surfCellFacesOnProcPatch[facei]];
            }

//          Pout<< "Sent at time = " << mesh_.time().value()
//              << ": surfCellFacesOnProcPatch = " << surfCellFacesOnProcPatch
//              << endl;
//          Pout<< "Sent at time = " << mesh_.time().value()
//              << ": dVfPatch = " << dVfPatch << endl;


//          toNbr << pFlux << surfCellFacesOnProcPatch << dVfPatch;

            toNbr << surfCellFacesOnProcPatch << dVfPatch;
        }

        pBufs.finishedSends();


        // Receive and combine
        forAll(procPatchLabels_, patchLabeli)
        {
            const label patchi = procPatchLabels_[patchLabeli];

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchi]);

            UIPstream fromNeighb(procPatch.neighbProcNo(), pBufs);
            DynamicList<label> faceIDs(100);
            DynamicList<scalar> nbrdVfs(100);
//          const scalarField nbrFlux(fromNeighb);

/*
            scalarField nbrFlux(procPatch.size());

            fromNeighb >> nbrFlux >> faceIDs >> nbrdVfs;
*/

            fromNeighb >> faceIDs >> nbrdVfs;

//          Pout<< "Received at time = " << mesh_.time().value()
//              << ": surfCellFacesOnProcPatch = " << faceIDs << endl;
//          Pout<< "Received at time = " << mesh_.time().value()
//              << ": dVfPatch = " << nbrdVfs << endl;

            // Combine fluxes
            scalarField& localFlux = dVf.boundaryFieldRef()[patchi];

/*
            const scalarField& phib = phi.boundaryField()[patchi];
            forAll(nbrFlux, facei)
            {
                if (phib[facei] < 0)
                {
                    localFlux[facei] = -nbrFlux[facei];
                }
            }
*/

            forAll(faceIDs, i)
            {
                const label facei = faceIDs[i];
                localFlux[facei] = - nbrdVfs[i];
                if (mag(localFlux[facei] + nbrdVfs[i]) > 10*SMALL)
                {
                    Pout<< "localFlux[facei] = " << localFlux[facei]
                        << " and nbrdVfs[i] = " << nbrdVfs[i]
                        << " for facei = " << facei << endl;
                }
            }
        }
/*
        // Write out results for checking
        forAll(procPatchLabels_, patchLabeli)
        {
            const label patchi = procPatchLabels_[patchLabeli];
            const scalarField& localFlux = dVf.boundaryField()[patchi];
            Pout<< "time = " << mesh_.time().value() << ": localFlux = "
                << localFlux << endl;
        }
*/
        // Reinitialising list used for minimal parallel communication
        forAll(surfaceCellFacesOnProcPatches_, patchi)
        {
            surfaceCellFacesOnProcPatches_[patchi].clear();
        }
    }
}


void Foam::isoAdvection::checkIfOnProcPatch(const label facei)
{
    if (!mesh_.isInternalFace(facei))
    {
        const polyBoundaryMesh& patches = mesh_.boundaryMesh();
        const label patchi = patches.whichPatch(facei);

        if
        (
            isA<processorPolyPatch>(patches[patchi])
         && patches[patchi].size() > 0
        )
        {
            const label patchFacei = patches[patchi].whichFace(facei);
            surfaceCellFacesOnProcPatches_[patchi].append(patchFacei);
//            Pout<< "patchFacei = " << patchFacei << " for facei = " << facei
//                << " on patchi = " << patchi << endl;
        }
    }
}


void Foam::isoAdvection::advect()
{

    // Interpolating alpha1 cell centre values to mesh points (vertices)
    ap_ = vpi_.interpolate(alpha1_.oldTime());

    // Initialising dVf with upwind values
    // i.e. phi[facei]*alpha1[upwindCell]*dt
    dVf_ =
        upwind<scalar>(mesh_, phi_).flux(alpha1_.oldTime())
       *mesh_.time().deltaT();

    // Do the isoAdvection on surface cells
    timeIntegratedFlux();

    // Synchronize processor patches
    syncProcPatches(dVf_, phi_);

    // Adjust dVf for unbounded cells
    limitFluxes();

    // Advect the free surface
    alpha1_ = alpha1_.oldTime() - fvc::surfaceIntegrate(dVf_);
    alpha1_.correctBoundaryConditions();
}


void Foam::isoAdvection::getSurfaceCells(cellSet& surfCells) const
{
    surfCells.clear();
    forAll(surfCells_, i)
    {
        surfCells.insert(surfCells_[i]);
    }
}


void Foam::isoAdvection::getBoundedCells(cellSet& boundCells) const
{
    boundCells.clear();
    forAll(cellIsBounded_, i)
    {
        if (cellIsBounded_[i])
        {
            boundCells.insert(i);
        }
    }
}


// ************************************************************************* //
