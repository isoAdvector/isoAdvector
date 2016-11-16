/*---------------------------------------------------------------------------*\

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
    along with IsoAdvector.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "isoAdvection.H"
#include "interpolationCellPoint.H"
#include "fvcSurfaceIntegrate.H"
#include "upwind.H"
#include "pointMesh.H"

// * * * * * * * * * * * * * * * *Optional debug * * * * * * * * * * * * * * //

#ifdef DETAILS2LOG
#define isoDebug(x) x
#else
#define isoDebug(x)
#endif


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
    alpha1In_(alpha1.oldTime().internalField()), // Need to reorganise
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
    nAlphaBounds_(dict_.lookupOrDefault<label>("nAlphaBounds", 1)),
    vof2IsoTol_(dict_.lookupOrDefault<scalar>("vof2IsoTol", 1e-8)),
    surfCellTol_(dict_.lookupOrDefault<scalar>("vof2IsoTol", 1e-8)),

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
        // communications get tangled)
        mesh_.C();

        // Get boundary mesh and resize the list for parallel comms
        const polyBoundaryMesh& patches = mesh_.boundaryMesh();
        surfaceCellFacesOnProcPatches_.resize(patches.size());

        // Append all processor patch labels to the list
        forAll(patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].size() > 0
            )
            {
                procPatchLabels_.append(patchI);
            }
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::isoAdvection::timeIntegratedFlux()
{
    isoDebug(Info << "Enter timeIntegratedFlux" << endl;)

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
    const scalarField& phiIn = phi_.internalField();
    const scalarField& magSfIn = mesh_.magSf().internalField();
    scalarField& dVfIn = dVf_.internalField();

    // Get necessary mesh data
    const labelListList& CP = mesh_.cellPoints();
    const labelListList& CC = mesh_.cellCells();
    const cellList& meshCells = mesh_.cells();
    const unallocLabelList& own = mesh_.owner();
    const unallocLabelList& nei = mesh_.neighbour();

    // Loop through cells
    forAll (alpha1In_, cellI)
    {
        if (isASurfaceCell(cellI))
        {
            // This is a surface cell, increment the counter, append and mark
            // the cell
            nSurfaceCells++;
            surfCells_.append(cellI);
            checkBounding_[cellI] = true;

            isoDebug
            (
                Info << "\n------------ Cell " << cellI << " with alpha1 = "
                    << alpha1In_[cellI] << " and 1-alpha1 = "
                    << 1.0-alpha1In_[cellI] << " ------------"
                    << endl;
            )

            // Calculate isoFace centre x0, normal n0 at time t
            label maxIter = 100; // NOTE: make it a debug switch

            // Calculate cell status (-1: cell is fully below the isosurface, 0:
            // cell is cut, 1: cell is fully above the isosurface)
            label cellStatus = isoCutCell_.vofCutCell
            (
                cellI,
                alpha1In_[cellI],
                vof2IsoTol_,
                maxIter
            );

//            Info << "1 - f0 = " << 1 - f0 << " for cell " << cellI << endl;

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
                    WarningIn
                    (
                        "void Foam::isoAdvection::timeIntegratedFlux()"
                    ) << "mag(n0) = " << mag(n0)
                      << " < 1e-6*minMagSf_ for cell " << cellI << endl;

                    // Initialise minimum and maximum values
                    scalar fMin = GREAT;
                    scalar fMax = -GREAT;

                    // Get cell points
                    const labelList& cellPts = CP[cellI];

                    // Calculate min and max values of the subset
                    subSetExtrema(ap_, cellPts, fMin, fMax);

                    scalar fInside  = 0;
                    if (alpha1In_[cellI] > 0.5) // Note: changed from >= to >
                    {
                        // Note: make a tolerance
                        fInside =  fMin + 1e-3*(fMax - fMin);
                    }
                    else
                    {
                        fInside =  fMax - 1e-3*(fMax - fMin);
                    }

//                    scalar fInside = f0 + sign(alpha1In_[cellI]-0.5)*1e-3;

                    // Calculate sub cell and initialise the normal with face
                    // area vector
                    isoCutCell_.calcSubCell(cellI, fInside);
                    n0 = isoCutCell_.isoFaceArea();
                }

                if (mag(n0) > 1e-6*minMagSf_)
                {
                    // Normalise the vector
                    isoDebug(Info << "Normalising n0: " << n0 << endl;)
                    n0 /= mag(n0);
                }
                else
                {
                    WarningIn
                    (
                        "void Foam::isoAdvection::timeIntegratedFlux()"
                    )   << "mag(n0) = " << mag(n0)
                        << " < 1e-6*minMagSf_ for cell " << cellI
                        << " with alpha1 = " << alpha1In_[cellI]
                        << ", 1-alpha1 = " << 1.0-alpha1In_[cellI]
                        << " and f0 = " << f0 << endl;

                    // Normalise the vector with stabilisation
                    n0 /= (mag(n0) + SMALL);
                    Info<< "After normalisation: mag(n0) = "
                        << mag(n0) << endl;
                }

                // Get the speed of the isoface by interpolating velocity and
                // dotting it with isoface normal
                const scalar Un0 = UInterp.interpolate(x0, cellI) & n0;

                isoDebug
                (
                    Info << "calcIsoFace gives initial surface: \nx0 = " << x0
                        << ", \nn0 = " << n0 << ", \nf0 = " << f0
                        << ", \nUn0 = " << Un0 << endl;
                )

                // Estimating time integrated water flux through each downwind
                // face
                const cell& cellFaces = meshCells[cellI];
                forAll (cellFaces, fi)
                {
                    // Get current face index
                    const label faceI = cellFaces[fi];

                    // Check if the face is internal face
                    if (mesh_.isInternalFace(faceI))
                    {
                        bool isDownwindFace = false;
                        label otherCell = -1;

                        // Check if the cell is owner
                        if (cellI == own[faceI])
                        {
                            // Note: writing phiIn > 0.0 instead of 10*SMALL
                            if (phiIn[faceI] > 0.0)
                            {
                                isDownwindFace = true;
                            }

                            // Other cell is neighbour
                            otherCell = nei[faceI];
                        }
                        else //Cell is the neighbour
                        {
                            // Note: writing phIn < 0.0 instead of -10*SMALL
                            if (phiIn[faceI] < 0.0)
                            {
                                isDownwindFace = true;
                            }

                            // Other cell is the owner
                            otherCell = own[faceI];
                        }

                        // Calculate time integrated flux if this is a downwind
                        // face
                        if (isDownwindFace)
                        {
//                            Info<< "Setting value for internal face " << faceI
//                              << endl;
                            dVfIn[faceI] = timeIntegratedFlux
                            (
                                faceI,
                                x0,
                                n0,
                                Un0,
                                f0,
                                dt,
                                phiIn[faceI],
                                magSfIn[faceI]
                            );
                        }

                        // We want to check bounding of neighbour cells to
                        // surface cells as well.
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
                        bsFaces_.append(faceI);
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
    const surfaceScalarField::GeometricBoundaryField& phib =
        phi_.boundaryField();
    const surfaceScalarField::GeometricBoundaryField& magSfb =
        mesh_.magSf().boundaryField();
    surfaceScalarField::GeometricBoundaryField& dVfb = dVf_.boundaryField();

    // Loop through boundary surface faces
    forAll(bsFaces_, fi)
    {
        // Get boundary face index (in the global list)
        const label fLabel = bsFaces_[fi];
//        Info << "Boundary face: " << fLabel << endl;

        // Get necesary mesh data
        const fvBoundaryMesh& boundaryMesh = mesh_.boundary();
        const polyBoundaryMesh& pBoundaryMesh = mesh_.boundaryMesh();

        // Get necessary labels
        // Note: consider optimisation since whichPatch is expensive
        const label patchI = pBoundaryMesh.whichPatch(fLabel);
        const label start = boundaryMesh[patchI].patch().start();
        const label size = boundaryMesh[patchI].size();

        if (size > 0)
        {
            // Get patch local label
            const label fLocal = fLabel - start;
            const scalar& phiP = phib[patchI][fLocal];

            if (phiP > 0) // Note: changed from phiP > 10*SMALL
            {
                const scalar& magSf = magSfb[patchI][fLocal];

                dVfb[patchI][fLocal] = timeIntegratedFlux
                (
                    fLabel,
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
                checkIfOnProcPatch(fLabel);
            }
        }
    }

    // Print out number of surface cells
    Info<< "Number of isoAdvector surface cells = "
        << returnReduce(nSurfaceCells, sumOp<label>()) << endl;
}


Foam::scalar Foam::isoAdvection::timeIntegratedFlux
(
    const label fLabel,
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

    isoDebug(Info << "Enter timeIntegratedFlux for face " << fLabel << endl;)

    // Treating rare cases where isoface normal is not calculated properly
    if (mag(n0) < 0.5)
    {
        scalar alphaf = 0.0;
        scalar waterInUpwindCell = 0.0;

        if (phi > 0 || !mesh_.isInternalFace(fLabel))
        {
            label upwindCell = mesh_.owner()[fLabel];
            alphaf = alpha1In_[upwindCell];
            waterInUpwindCell = alphaf*mesh_.V()[upwindCell];
        }
        else
        {
            label upwindCell = mesh_.neighbour()[fLabel];
            alphaf = alpha1In_[upwindCell];
            waterInUpwindCell = alphaf*mesh_.V()[upwindCell];
        }

        isoDebug
        (
            WarningIn
            (
                "Foam::scalar Foam::isoAdvection::timeIntegratedFlux(...)"
            ) << "mag(n0) = " << mag(n0)
              << " so timeIntegratedFlux calculates dVf from upwind"
              << " cell alpha value: " << alphaf << endl;
        )

        return min(alphaf*phi*dt, waterInUpwindCell);
    }

    // Find sorted list of times where the isoFace will arrive at face points
    // given initial position x0 and velocity Un0*n0

    // Get points for this face
    const face& pLabels = mesh_.faces()[fLabel];

    // Note: changed to direct access to points from the face
    const pointField fPts(pLabels.points(mesh_.points()));
    const label nPoints = fPts.size();
    
    scalarField pTimes(fPts.size());
    if (mag(Un0) > 1e-12) // Note: tolerances
    {
        // Here we estimate time of arrival to the face points from their normal
        // distance to the initial surface and the surface normal velocity

        pTimes = ((fPts - x0) & n0)/Un0;
        
        scalar dVf = 0;
        
        //Check if pTimes changes direction more than twice when looping face
        label nShifts = 0;
        forAll(pTimes, pi)
        {
            label oldEdgeSign = sign(pTimes[(pi + 1) % nPoints] - pTimes[pi]);
            label newEdgeSign = sign(pTimes[(pi + 2) % nPoints] - pTimes[(pi + 1) % nPoints]);
            if (newEdgeSign != oldEdgeSign)
            {
                nShifts++;
            }
        }
        if (nShifts == 2)
        {
            dVf = phi/magSf*timeIntegratedArea(fPts,pTimes,dt,magSf,Un0);
        }
        else if (nShifts > 2) //triangle decompose the face
        {
            pointField fPts_tri(3);
            scalarField pTimes_tri(3);
            fPts_tri[0] = mesh_.Cf()[fLabel];
            pTimes_tri[0] = ((fPts_tri[0] - x0) & n0)/Un0;
            for (label pi = 0; pi < nPoints; pi++)
            {
                fPts_tri[1] = fPts[pi];
                pTimes_tri[1] = pTimes[pi];
                fPts_tri[2] = fPts[(pi + 1) % nPoints];
                pTimes_tri[2] = pTimes[(pi + 1) % nPoints];
                const scalar magSf_tri = mag(0.5*(fPts_tri[2] - fPts_tri[0]) ^ (fPts_tri[1] - fPts_tri[0]));
                const scalar phi_tri = phi*magSf_tri/magSf;
                dVf += phi_tri/magSf_tri*timeIntegratedArea(fPts_tri,pTimes_tri,dt,magSf_tri,Un0);
            }
//            WarningIn
//            (
//                "Foam::scalar Foam::isoAdvection::timeIntegratedFlux(...)"
//            )   << "Warning: nShifts = " << nShifts << " on face " << fLabel 
//                << " owned by cell " << mesh_.owner()[fLabel] << endl;            
        }
        else
        {
            WarningIn
            (
                "Foam::scalar Foam::isoAdvection::timeIntegratedFlux(...)"
            )   << "Warning: nShifts = " << nShifts << " on face " << fLabel
                << " with pTimes = " << pTimes << " owned by cell " 
                << mesh_.owner()[fLabel] << endl;
        }

        return dVf;
    }
    else
    {
        // Un0 is almost zero and isoFace is treated as stationary
        isoCutFace_.calcSubFace(fLabel, f0);
        scalar alphaf = mag(isoCutFace_.subFaceArea()/magSf);
        isoDebug
        (
            Info                << endl;
            WarningIn
            (
                "Foam::scalar Foam::isoAdvection::timeIntegratedFlux(...)"
            ) << "Warning: Un0 is almost zero (" << Un0
              << ") - calculating dVf on face " << fLabel
              << " using subFaceFraction giving alphaf = " << alphaf
              << endl;
        )

        return phi*dt*alphaf;
    }
}


Foam::scalar Foam::isoAdvection::timeIntegratedArea
(
    const pointField& fPts,
    const scalarField& pTimes,
    const scalar dt,
    const scalar magSf,
    const scalar Un0
)
{
    isoDebug(Info << "Enter timeIntegratedArea with fPts = " << fPts << endl;)

    // Initialise time integrated area
    scalar tIntArea = 0.0;

    // Sorting face vertex encounter time list
    scalarList sortedTimes(pTimes);
    sort(sortedTimes);

    isoDebug(Info << "sortedTimes = " << sortedTimes << endl;)

    // Dealing with case where face is not cut by surface during time interval
    // [0, dt] because face was already passed by surface
    if (sortedTimes.last() < SMALL) // Note: changed from <= 0.0
    {
        // All cuttings in the past
        isoDebug(Info << "All cuttings in the past" << endl;)

        // If all face cuttings were in the past and cell is filling up (Un0>0)
        // then face must be full during whole time interval
        tIntArea = magSf*dt*pos(Un0);
        return tIntArea;
    }

    // Dealing with case where face is not cut by surface during time interval
    // [0, dt] because dt is too small for surface to reach closest face point
    if (sortedTimes.first() > dt) // Note: changed from >= dt
    {
        // All cuttings in the future
        isoDebug(Info << "All cuttings in the future" << endl;)

        // If all cuttings are in the future but non of them within [0, dt] then
        // if cell is filling up (Un0 > 0) face must be empty during whole time
        // interval
        tIntArea = magSf*dt*(1-pos(Un0));
        return tIntArea;
    }

    // Cutting sortedTimes at 0 and dt and sorting out duplicates
    DynamicScalarList t(sortedTimes.size() + 2);
    t.append(0);
    scalar smallTime = max(1e-6*dt, 10*SMALL); // Note: tolerances
    forAll(sortedTimes, ti)
    {
        const scalar& curTime = sortedTimes[ti];

        if
        (
            smallTime < curTime
         && curTime < dt - smallTime
         && mag(curTime - t.last()) > smallTime
        )
        {
            t.append(curTime);
        }
    }

    // Finally append time step
    t.append(dt);

    isoDebug(Info << "Cutting sortedTimes at 0 and dt: t = " << t << endl;)

//    Info << "times = " << t << endl;

    // Get information whether the face is uncut during first and last interval
    bool faceUncutInFirstInterval(sortedTimes.first() > 0.0);
    bool faceUncutInLastInterval(sortedTimes.last() < dt);

    // Dealing with cases where face is cut at least once during time interval
    // [0, dt]
    label nt = 0; // Sub time interval counter
    scalar initialArea = 0.0; //Submerged area at time

    // Special treatment for first time interval if face is uncut during this
    if (faceUncutInFirstInterval)
    {
        // If Un0 > 0 cell is filling up - hence if face is cut at a later time
        // but not initially it must be initially empty
        tIntArea = magSf*(t[nt + 1] - t[nt])*(1.0 - pos(Un0));
        initialArea = magSf*(1.0 - pos(Un0));
        isoDebug
        (
            Info<< "faceUncutInFirstInterval, so special treatment for"
                << " first time interval: [" << t[nt] << ", " << t[nt+1]
                << "] giving tIntArea = " << tIntArea << endl;
        )
        nt++;
    }
    else // Calculate initialArea if face is initially cut
    {
        isoDebug
        (
            Info<< "face is initially cut, so finding initial area, pTimes = "
                << pTimes << ", Un0 = " << Un0 << endl;
        )

        isoCutFace_.calcSubFace(fPts, -sign(Un0)*pTimes, 0.0);
        initialArea = mag(isoCutFace_.subFaceArea());
    }

    isoDebug
    (
        Info<< "InitialArea for next time step corresponds to face phase"
            << " fraction a0 = " << initialArea/magSf << " where |Sf| = "
            << magSf << " was used." << endl;
    )
    while (nt < t.size() - (1 + faceUncutInLastInterval))
    {
        // Note: please think of grouping some of the lines in logical fashion
        DynamicPointList cutPoints1(3), cutPoints2(3);
        isoCutFace_.cutPoints(fPts, pTimes, t[nt], cutPoints1);
        isoCutFace_.cutPoints(fPts, pTimes, t[nt + 1], cutPoints2);
        scalar quadArea, intQuadArea;
        quadAreaCoeffs(cutPoints1, cutPoints2, quadArea, intQuadArea);
        // Integration of area(t) = A*t^2+B*t from t = 0 to 1.
        scalar integratedQuadArea = sign(Un0)*intQuadArea;
        tIntArea += (t[nt + 1] - t[nt])*(initialArea + integratedQuadArea);
        // Adding quad area to submerged area
        initialArea += sign(Un0)*quadArea;

        isoDebug
        (
            Info<< "Integrating area for " << nt + 1 << "'th time interval: ["
                << t[nt] << ", " << t[nt + 1] << "] giving tIntArea = "
                << tIntArea << " and a0 = " << initialArea/magSf << endl;
        )
        nt++;
    }

    //Special treatment for last time interval if face is uncut during this
    if (faceUncutInLastInterval && (nt + 1) < t.size())
    {
        isoDebug
        (
            Info<< "faceUncutInLastInterval, so special treatment for last ("
                << nt+1 << "'th) time interval: [" << t[nt] << ", " << t[nt+1]
                << "]" << endl;
        )
        // If face is cut at some intermediate time but not at last time, then
        // if Un0 > 0 (cell filling up) face must be filled at last time
        // interval.
        tIntArea += magSf*(t[nt + 1] - t[nt])*pos(Un0);
    }

    return tIntArea;
}


void Foam::isoAdvection::getDownwindFaces
(
    const label cellI,
    DynamicLabelList& downwindFaces
)
{
    // Note: this function is called within a loop, consider passing mesh owners
    // as arguments

    isoDebug(Info << "Enter getDownwindFaces " << endl;)

    // Get necessary mesh data and cell information
    const unallocLabelList& own = mesh_.owner();
    const cellList& cells = mesh_.cells();
    const cell& fLabels = cells[cellI];

    // Check all faces of the cell
    forAll(fLabels, fi)
    {
        // Get face and corresponding flux
        const label fLabel = fLabels[fi];
        const scalar& phi = faceValue(phi_, fLabel);

        if (own[fLabel] == cellI)
        {
            if (phi > 0.0) // Note: changed from > 10*SMALL
            {
                downwindFaces.append(fLabel);
            }
        }
        else if (phi < 0.0) // Note: changed from < -10*SMALL
        {
            downwindFaces.append(fLabel);
        }
    }

    downwindFaces.shrink();
}


void Foam::isoAdvection::quadAreaCoeffs
(
    const DynamicPointList& pf0,
    const DynamicPointList& pf1,
    scalar& quadArea,
    scalar& intQuadArea
)
{
    isoDebug(Info << "Enter quadAreaCoeffs" << endl;)
//    isoDebug
//    (
//        Info<< "Face " << fLabel << " was cut at " << pf0 << " by f0 = "
//            << f0 << " and at " << pf1 << " by " << " f1 = " << f1
//            << endl;
//    )

    const label np0 = pf0.size();
    const label np1 = pf1.size();

    scalar alpha = 0.0;
    scalar beta = 0.0;
    quadArea = 0.0;
    intQuadArea = 0.0;

    if (np0 > 0 && np1 > 0)
    {
        // Defining quadrilateral vertices
        vector A(pf0[0]);
        vector C(pf1[0]);
        vector B(vector::zero);
        vector D(vector::zero);

        // Triangle cases
        if (np0 == 2 && mag(pf0[0] - pf0[1]) > SMALL)
        {
            B = pf0[1];
        }
        else
        {
            // Note: tolerances
            B = A + 1e-4*(pf1[1] - pf1[0]);
            if (np0  != 1)
            {
                WarningIn
                (
                    "void Foam::isoAdvection::quadAreaCoeffs(...)"
                ) << "Vertex face was cut at pf0 = " << pf0 << endl;
            }
        }

        if (np1 == 2 && mag(pf1[0] - pf1[1]) > SMALL)
        {
            D = pf1[1];
        }
        else
        {
            // Note: tolerances
            D = C + 1e-4*(A - B);
            if (np1 != 1)
            {
                WarningIn
                (
                    "void Foam::isoAdvection::quadAreaCoeffs(...)"
                ) << "Vertex face was cut at pf1 = " << pf1 << endl;
            }
        }

        // Defining local coordinates for area integral calculation
        vector xhat = B - A;
        xhat /= mag(xhat);
        vector yhat = (D-A);
        yhat -= ((yhat & xhat) * xhat);
        yhat /= mag(yhat);

//        isoDebug
//        (
//            Info<< "xhat = " << xhat << ", yhat = " << yhat << ", zhat = "
//                << zhat << ". x.x = " << (xhat & xhat) << ", y.y = "
//                << (yhat & yhat) <<", z.z = " << (zhat & zhat) << ", x.y = "
//                << (xhat & yhat) << ", x.z = " << (xhat & zhat) << ", y.z = "
//                << (yhat & zhat) << endl;
//         )


        // Swapping pf1 points if pf0 and pf1 point in same general direction
        // (because we want a quadrilateral ABCD where pf0 = AB and pf1 = CD)
        if (((B - A) & (D - C)) > 0)
        {
//            Info << "Swapping C and D" << endl;
            vector tmp = D;
            D = C;
            C = tmp;
        }

//        Info<< "A = " << A << ", B = " << B << ", C = " << C << ", D = "
//            << D << endl;
        scalar Bx = mag(B - A);
        scalar Cx = (C - A) & xhat;
        scalar Cy = mag((C - A) & yhat);
        scalar Dx = (D - A) & xhat;
        scalar Dy = mag((D - A) & yhat);

//      area = ((Cx-Bx)*Dy-Dx*Cy)/6.0 + 0.25*Bx*(Dy+Cy);
        alpha = 0.5*((Cx - Bx)*Dy - Dx*Cy);
        beta = 0.5*Bx*(Dy + Cy);
        quadArea = alpha + beta;

        // Possible bug: 1/3 probably evaluated as 0?
//        intQuadArea = (1/3)*alpha + 0.5*beta;
        intQuadArea = alpha/3.0 + 0.5*beta;

//         Info<< "Bx = " << Bx << ", Cx = " << Cx << ", Cy = " << Cy
//             << ", Dx = " << Dx << ", Dy = " << Dy << ", alpha = " << alpha
//             << ", beta = " << beta << endl;

        //area(t) = A*t^2+B*t
        //integratedArea = A/3+B/2
    }
    else
    {
        Info<< "Vertex face was cut at " << pf0
            << " and at " << pf1 << " by " << endl;
    }
}


void Foam::isoAdvection::subSetExtrema
(
    const scalarField& f,
    const labelList& labels,
    scalar& fMin,
    scalar& fMax
)
{
    fMin = VGREAT;
    fMax = -VGREAT;

    forAll(labels,pi)
    {
        scalar fp = f[labels[pi]];

        if (fp < fMin)
        {
            fMin = fp;
        }

        if (fp > fMax)
        {
            fMax = fp;
        }
    }
}


void Foam::isoAdvection::limitFluxes()
{
    // Get time step size
    const scalar dt = mesh_.time().deltaT().value();

//    scalarField alphaNew = alpha1In_ - fvc::surfaceIntegrate(dVf);
    scalar aTol = 1.0e-12; // Note: tolerances
    scalar maxAlphaMinus1 = 1; // max(alphaNew - 1);
    scalar minAlpha = -1; // min(alphaNew);
    label nUndershoots = 20; // sum(neg(alphaNew + aTol));
    label nOvershoots = 20; // sum(pos(alphaNew - 1 - aTol));
    cellIsBounded_ = false;

    // Loop number of bounding steps
    for (label n = 0; n < nAlphaBounds_; n++)
    {
        Info<< "Running bounding number " << n + 1 << " of time "
            << mesh_.time().value() << endl;

        if (maxAlphaMinus1 > aTol) // Note: tolerances
        {
            isoDebug(Info << "Bound from above... " << endl;)

//          scalarField dVfcorrected = dVf.internalField();

            surfaceScalarField dVfcorrected("dVfcorrected", dVf_);
            DynamicLabelList correctedFaces(3*nOvershoots);
            boundFromAbove(alpha1In_, dVfcorrected, correctedFaces);

            forAll(correctedFaces, fi)
            {
                label fLabel = correctedFaces[fi];

                // Change to treat boundaries consistently
                faceValue(dVf_, fLabel, faceValue(dVfcorrected, fLabel));
            }

            syncProcPatches(dVf_, phi_);
        }

        if (minAlpha < -aTol) // Note: tolerances
        {
            isoDebug(Info << "Bound from below... " << endl;)

            scalarField alpha2 = 1.0 - alpha1In_;
            surfaceScalarField dVfcorrected
            (
                "dVfcorrected",
                phi_*dimensionedScalar("dt", dimTime, dt) - dVf_
            );
//            dVfcorrected -= dVf; // phi_ and dVf have same sign and dVf is the
//                                 // portion of phi_*dt that is water.
            // If phi_ > 0 then dVf > 0 and mag(phi_*dt-dVf) < mag(phi_*dt) as
            // it should.
            // If phi_ < 0 then dVf < 0 and mag(phi_*dt-dVf) < mag(phi_*dt) as
            // it should.
            DynamicLabelList correctedFaces(3*nUndershoots);
            boundFromAbove(alpha2, dVfcorrected, correctedFaces);

            forAll(correctedFaces,fi)
            {
                label fLabel = correctedFaces[fi];

                //Change to treat boundaries consistently
                scalar phi = faceValue(phi_, fLabel);
                scalar dVcorr = faceValue(dVfcorrected, fLabel);
                faceValue(dVf_, fLabel, phi*dt - dVcorr);
            }

            syncProcPatches(dVf_, phi_);
        }

//        //Check if still unbounded
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
    DynamicLabelList& correctedFaces
)
{
    // Get time step size
    const scalar dt = mesh_.time().deltaT().value();

    correctedFaces.clear();
    scalar aTol = 10*SMALL; // Note: tolerances

    // Get necessary mesh data
    const scalarField& meshV = mesh_.V();
    const cellList& meshCells = mesh_.cells();

    // Loop through alpha cell centred field
    forAll(alpha1, cellI)
    {
        if (checkBounding_[cellI])
        {
            const scalar& Vi = meshV[cellI];
            scalar alpha1New = alpha1[cellI] - netFlux(dVf, cellI)/Vi;
            scalar alphaOvershoot = alpha1New - 1.0;
            scalar fluidToPassOn = alphaOvershoot*Vi;
            label nFacesToPassFluidThrough = 1;

            bool firstLoop = true;
            // First try to pass surplus fluid on to neighbour cells that are
            // not filled and to which dVf < phi*dt
            while (alphaOvershoot > aTol && nFacesToPassFluidThrough > 0)
            {
                cellIsBounded_[cellI] = true;
                isoDebug
                (
                    Info<< "\n\nBounding cell " << cellI
                        << " with alpha overshooting " << alphaOvershoot
                        << endl;
                )

                // First find potential neighbour cells to pass surplus water to
                DynamicLabelList downwindFaces(meshCells[cellI].size());
                getDownwindFaces(cellI, downwindFaces);

                DynamicLabelList facesToPassFluidThrough(downwindFaces.size());
                DynamicScalarList dVfmax(downwindFaces.size());
                DynamicScalarList phi(downwindFaces.size());

                scalar dVftot = 0.0;
                nFacesToPassFluidThrough = 0;

                isoDebug(Info << "downwindFaces: " << downwindFaces << endl;)

                forAll(downwindFaces, fi)
                {
                    const label fLabel = downwindFaces[fi];
                    scalar maxExtraFaceFluidTrans = 0.0;
                    const scalar phif = faceValue(phi_, fLabel);
                    const scalar dVff = faceValue(dVf, fLabel);
                    maxExtraFaceFluidTrans = mag(phif*dt - dVff);

                    // dVf has same sign as phi and so if phi > 0 we have
                    // mag(phi_[fLabel]*dt) - mag(dVf[fLabel]) = phi_[fLabel]*dt
                    // - dVf[fLabel]
                    // If phi < 0 we have mag(phi_[fLabel]*dt) -
                    // mag(dVf[fLabel]) = -phi_[fLabel]*dt - (-dVf[fLabel]) > 0
                    // since mag(dVf) < phi*dt
                    isoDebug
                    (
                        Info<< "downwindFace " << fLabel
                            << " has maxExtraFaceFluidTrans = "
                            << maxExtraFaceFluidTrans << endl;
                    )

                    if (maxExtraFaceFluidTrans/Vi > aTol) // Note: tolerances
//                    if ( maxExtraFaceFluidTrans/Vi > aTol &&
//                    mag(dVfIn[fLabel])/Vi > aTol ) //Last condition may be
//                    important because without this we will flux through uncut
//                    downwind faces
                    {
                        facesToPassFluidThrough.append(fLabel);
                        phi.append(phif);
                        dVfmax.append(maxExtraFaceFluidTrans);
                        dVftot += mag(phif*dt);
                    }
                }

                isoDebug
                (
                    Info<< "\nfacesToPassFluidThrough: "
                        << facesToPassFluidThrough << ", dVftot = "
                        << dVftot << " m3 corresponding to dalpha = "
                        << dVftot/Vi << endl;
                )

                forAll(facesToPassFluidThrough, fi)
                {
                    const label fLabel = facesToPassFluidThrough[fi];
                    scalar fluidToPassThroughFace = fluidToPassOn*
                        mag(phi[fi]*dt)/dVftot;

                    nFacesToPassFluidThrough +=
                        pos(dVfmax[fi] - fluidToPassThroughFace);

                    fluidToPassThroughFace =
                        min(fluidToPassThroughFace, dVfmax[fi]);

                    scalar dVff = faceValue(dVf, fLabel);
                    dVff += sign(phi[fi])*fluidToPassThroughFace;
                    faceValue(dVf, fLabel, dVff);

                    if (firstLoop)
                    {
                        checkIfOnProcPatch(fLabel);
                        correctedFaces.append(fLabel);
                    }
                }

                firstLoop = false;
                alpha1New = alpha1[cellI] - netFlux(dVf, cellI)/Vi;
                alphaOvershoot = alpha1New - 1.0;
                fluidToPassOn = alphaOvershoot*Vi;

                isoDebug
                (
                    Info<< "\nNew alpha for cell " << cellI << ": "
                        << alpha1New << endl;
                )
            }
        }
    }

    isoDebug(Info << "correctedFaces = " << correctedFaces << endl;)
}


Foam::scalar Foam::isoAdvection::netFlux
(
    const surfaceScalarField& dVf,
    const label cellI
)
{
    scalar dV = 0.0;

    // Get face label
    const labelList& fLabels = mesh_.cells()[cellI];

    // Get mesh data
    const unallocLabelList& own = mesh_.owner();

    forAll (fLabels, fi)
    {
        const label fLabel = fLabels[fi];
        const scalar dVff = faceValue(dVf, fLabel);

        if (own[fLabel] == cellI)
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
    // Get list of patches
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    if (Pstream::parRun())
    {
        // Send data
        forAll(procPatchLabels_, patchLabelI)
        {
            // Get current patch label
            const label patchI = procPatchLabels_[patchLabelI];

            // Get reference to current processor patch
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchI]);

            // Get flux for this patch
            const scalarField& pFlux = dVf.boundaryField()[patchI];

            // Get the list of current surfaceCell faces on this processor patch
            const labelList& surfCellFacesOnProcPatch =
                surfaceCellFacesOnProcPatches_[patchI];

            // Calculate the field that will be sent to the other side
            scalarList dVfPatch(surfCellFacesOnProcPatch.size());
            forAll(dVfPatch, i)
            {
                dVfPatch[i] = pFlux[surfCellFacesOnProcPatch[i]];
            }

            // Send data to neighbouring processor
            OPstream toNbr(Pstream::blocking, procPatch.neighbProcNo());
            toNbr << surfCellFacesOnProcPatch << dVfPatch;
        }

        // Receive data and combine
        forAll(procPatchLabels_, patchLabelI)
        {
            // Get current patch label
            const label patchI = procPatchLabels_[patchLabelI];

            // Get reference to current processor patch
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchI]);

            // Receive data
            labelList fLabels;
            scalarList nbrdVfs;
            IPstream fromNbr(Pstream::blocking, procPatch.neighbProcNo());
            fromNbr >> fLabels >> nbrdVfs;

            // Combine fluxes
            scalarField& localFlux = dVf.boundaryField()[patchI];

            forAll (fLabels, faceI)
            {
                // Get label of the face
                const label& fLabel = fLabels[faceI];
                localFlux[fLabel] = - nbrdVfs[faceI];

                if (mag(localFlux[fLabel] + nbrdVfs[faceI]) > 10*SMALL)
                {
                    Pout<< "localFlux[fLabel] = " << localFlux[fLabel]
                        << " and nbrdVfs[faceI] = " << nbrdVfs[faceI]
                        << " for fLabel = " << fLabel << endl;
                }
            }
        }

        // Reinitialising list used for minimal parallel communication
        forAll(surfaceCellFacesOnProcPatches_,patchI)
        {
            surfaceCellFacesOnProcPatches_[patchI].clear();
        }
    }
}


void Foam::isoAdvection::checkIfOnProcPatch
(
    const label faceI
)
{
    if (!mesh_.isInternalFace(faceI))
    {
        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        const label patchI = patches.whichPatch(faceI);

        if
        (
            isA<processorPolyPatch>(patches[patchI])
         && patches[patchI].size() > 0
        )
        {
            const label fLabel = patches[patchI].whichFace(faceI);
            surfaceCellFacesOnProcPatches_[patchI].append(fLabel);
//            Pout<< "fLabel = " << fLabel << " for faceI = " << faceI
//                << " on patchI = " << patchI << endl;
        }
    }
}


void Foam::isoAdvection::advect()
{
    isoDebug(Info << "Enter advect" << endl;)

    // Interpolating alpha1 cell centre values to mesh points (vertices)
    ap_ = vpi_.interpolate(alpha1_.oldTime());

    // Initialising dVf with upwind values, i.e.
    // phi[fLabel]*alpha1[upwindCell]*dt
    dVf_ = upwind<scalar>(mesh_, phi_).flux(alpha1_.oldTime())*
        mesh_.time().deltaT();

    // Do the isoAdvection on surface cells
    timeIntegratedFlux();

    // Syncronize processor patches
    syncProcPatches(dVf_, phi_);

    // Adjust dVf for unbounded cells
    limitFluxes();

    // Advect the free surface
    alpha1_ = alpha1_.oldTime() - fvc::surfaceIntegrate(dVf_);
    alpha1_.correctBoundaryConditions();
}


void Foam::isoAdvection::getSurfaceCells
(
    cellSet& surfCells
)
{
    surfCells.clear();
    forAll(surfCells_, i)
    {
        surfCells.insert(surfCells_[i]);
    }
}


void Foam::isoAdvection::getBoundedCells
(
    cellSet& boundCells
)
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
