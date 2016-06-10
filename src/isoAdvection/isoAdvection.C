/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "isoAdvection.H"
#include "volPointInterpolation.H"
#include "interpolationCellPoint.H"
#include "fvcSurfaceIntegrate.H"
#include "upwind.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

#ifdef DETAILS2LOG
#define isoDebug(x) x
#else
#define isoDebug(x)
#endif


Foam::isoAdvection::isoAdvection
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U
)
:
    mesh_(alpha1.mesh()),
    alpha1_(alpha1),
    alpha1In_(alpha1.internalField()),
    vpi_(mesh_),
    ap_(mesh_.nPoints(),0.0),
    phi_(phi),
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
		dimensionedScalar("vol", dimVol, 0)
    ),
    U_(U),
    nAlphaBounds_(1),
    vof2IsoTol_(1e-8),
    surfCellTol_(1e-8),
    isoCutCell_(mesh_,ap_),
    isoCutFace_(mesh_,ap_),
    procPatchLabels_(0),
    surfaceCellFacesOnProcPatches_(0),
    surfCells_(max(10, ceil(0.2*mesh_.nCells()))),
    cellIsBounded_(mesh_.nCells(),false),
    checkBounding_(mesh_.nCells(),false),
    bsFaces_(max(10, ceil(0.2*(mesh_.nFaces()-mesh_.nInternalFaces())))),
    bsx0_(bsFaces_.size()),
    bsn0_(bsFaces_.size()),
    bsUn0_(bsFaces_.size()),
    bsf0_(bsFaces_.size()),
    minMagSf_(gMin(mesh_.magSf()))
{
    
    //Reading isoAdvector controls from fvSolution if present
    const dictionary& isoAdvectorDict = mesh_.solutionDict().subDict("isoAdvector");
    nAlphaBounds_ = isoAdvectorDict.lookupOrDefault<label>("nAlphaBounds", 1);
    vof2IsoTol_ = isoAdvectorDict.lookupOrDefault<scalar>("vof2IsoTol", 1e-8);
    surfCellTol_ = isoAdvectorDict.lookupOrDefault<scalar>("surfCellTol", 1e-8);


    //Prepare lists used in parallel runs
    if (Pstream::parRun())
    {
        //Force calculation of cell centres and volumes (else parallel 
        //communication may crash)
        mesh_.C();

        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        surfaceCellFacesOnProcPatches_.resize(patches.size());

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


void Foam::isoAdvection::timeIntegratedFlux
(
    const scalar dt,
    surfaceScalarField& dVf
)
{
    isoDebug(Info << "Enter timeIntegratedFlux" << endl;)

    //Calculating weights for interpolation of U to the isoface face centre
    interpolationCellPoint<vector> UInterp(U_);
    
    isoDebug(Info << "Running through all surface cells" << endl;)

    //For each downwind face of each surface cell we "isoadvect" to find dVf
    label nSurfaceCells = 0;
    surfCells_.clear();
    checkBounding_ = false;
    bsFaces_.clear();
    bsx0_.clear();
    bsn0_.clear();
    bsUn0_.clear();
    bsf0_.clear();
    const scalarField& phiIn = phi_.internalField();
    const scalarField& magSfIn = mesh_.magSf().internalField();
    scalarField& dVfIn = dVf.internalField();
    
    forAll (alpha1In_, cellI)
    {
        if ( isASurfaceCell(cellI) )
        {
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

            //Calculate isoFace centre x0, normal n0 at time t
            label maxIter(100);
//            scalar f0 = isoCutCell_.vofCutCell(cellI, alpha1In_[cellI], vof2IsoTol_, maxIter);
            label cellStatus = isoCutCell_.vofCutCell2(cellI, alpha1In_[cellI], vof2IsoTol_, maxIter);
//            Info << "1 - f0 = " << 1 - f0 << " for cell " << cellI << endl;
            if (cellStatus == 0)
            {
                scalar f0 = isoCutCell_.isoValue();
                point x0 = isoCutCell_.isoFaceCentre();
                vector n0 = isoCutCell_.isoFaceArea();
                            
                //If cell almost full or empty isoFace may be undefined. 
                //Calculating normal by going a little into the cell.
                if ( mag(n0) < 1e-6*minMagSf_ ) 
                {
                    Info << "Warning: mag(n0) = " << mag(n0) 
                        << " < 1e-6*minMagSf_ for cell " << cellI << endl;
                    scalar fMin(GREAT), fMax(-GREAT);
                    const labelList& cellPts = mesh_.cellPoints()[cellI];
                    subSetExtrema(ap_, cellPts, fMin, fMax);
                    scalar fInside  = 0;
                    if (alpha1In_[cellI] >= 0.5)
                    {
                        fInside =  fMin + 1e-3*(fMax-fMin);
                    }
                    else 
                    {
                        fInside =  fMax - 1e-3*(fMax-fMin);                        
                    }
//                    scalar fInside = f0 + sign(alpha1In_[cellI]-0.5)*1e-3;
                    isoCutCell_.calcSubCell(cellI,fInside);
                    n0 = isoCutCell_.isoFaceArea();
                }    
                
                if ( mag(n0) > 1e-6*minMagSf_ )
                {
                    isoDebug(Info << "Normalising n0: " << n0 << endl;)
                    n0 /= mag(n0);
                }
                else
                {
                    Info << "Warning: mag(n0) = " << mag(n0)
                    << " < 1e-6*minMagSf_ for cell " << cellI << " with alpha1 = " 
                        << alpha1In_[cellI] << ", 1-alpha1 = " 
                        << 1.0-alpha1In_[cellI] << " and f0 = " << f0 << endl;
                    n0 /= (mag(n0) + SMALL);
                    Info << "After normalisation: mag(n0) = " << mag(n0) << endl;
                }

                //Interpolate velocity to isoFace centre
                vector U0 = UInterp.interpolate(x0,cellI);
                scalar Un0 = (U0 & n0);

                isoDebug
                (
                    Info << "calcIsoFace gives initial surface: \nx0 = " << x0 
                        << ", \nn0 = " << n0 << ", \nf0 = " << f0 << ", \nUn0 = " 
                        << Un0 << endl;
                )

                //Estimating time integrated water flux through each downwind face
                const labelList& cellFaces = mesh_.cells()[cellI];
                forAll (cellFaces, fi)
                {
                    const label faceI = cellFaces[fi];
                    if (mesh_.isInternalFace(faceI))
                    {
                        bool isDownwindFace = false;
                        label otherCell = -1;
                        if (cellI == mesh_.owner()[faceI])
                        {
    //                        if (phiIn[faceI] > 1e-12)
                            if (phiIn[faceI] > 10*SMALL)
                            {
                                isDownwindFace = true;
                            }
                            otherCell = mesh_.neighbour()[faceI];
                        }
                        else //cellI must be neighbour
                        {
    //                        if (phiIn[faceI] < -1e-12)
                            if (phiIn[faceI] < -10*SMALL)
                            {
                                isDownwindFace = true;
                            }
                            otherCell = mesh_.owner()[faceI];                            
                        }

                        if (isDownwindFace)
                        {
    //                        Info << "Setting value for internal face " << faceI << endl;
                            dVfIn[faceI] = timeIntegratedFlux(faceI, x0, n0, Un0,
                                f0, dt, phiIn[faceI], magSfIn[faceI]);
                        }
                        //We want to check bounding of neighbour cells to surface 
                        //cells as well:
                        checkBounding_[otherCell] = true;
                        //Also check neighbour neighbours
                        const labelList& nNeighbourCells = mesh_.cellCells()[otherCell];
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
                        checkIfOnProcPatch(faceI);
                    }
                }
            }
        }
    }
    
    //Setting dVf for boundary faces of surface cells
    if(bsFaces_.size() > 0)
    {
        label patchI = -1;
        label start = -1;
        label size = -1;

        forAll(bsFaces_, fi)
        {
            const label fLabel = bsFaces_[fi];
//            Info << "Boundary face: " << fLabel << endl;
            if (fLabel >= start + size)
            {
                patchI = mesh_.boundaryMesh().whichPatch(fLabel);
                start = mesh_.boundary()[patchI].start();
                size = mesh_.boundary()[patchI].size();
//                Info << "Changing to patch " << patchI << " with start = " 
//                    << start << " and size = " << size << endl;
            }
            if (size > 0)
            {
                const scalar phiP = phi_.boundaryField()[patchI][fLabel - start];
//                if (phiP > 1e-12)
                if (phiP > 10*SMALL)
                {
                    const scalar magSf = mesh_.magSf().boundaryField()[patchI][fLabel - start];
                    scalar& dVfP = dVf.boundaryField()[patchI][fLabel - start];
                    dVfP = timeIntegratedFlux(fLabel, bsx0_[fi],
                        bsn0_[fi], bsUn0_[fi], bsf0_[fi], dt,
                        phiP, magSf);
                }
            }
        }
    }
    Pout << "nSurfaceCells = " << nSurfaceCells << " out of " << alpha1In_.size() 
        << " cells" << endl;
}


bool Foam::isoAdvection::isASurfaceCell
(
    const label cellI
)
{
    return 
    (
        surfCellTol_ < alpha1In_[cellI] 
     && alpha1In_[cellI] < 1 - surfCellTol_ 
    );
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
    isoDebug(Info << "Enter timeIntegratedFlux for face " << fLabel << endl;)

    scalar dVf = 0; //Volume flowing through face in time interval [0,dt] to be calculated below

    //Treating rare cases where isoface normal is not calculated properly
    if (mag(n0) < .5)
    {
        scalar alphaf = 0.0;
        scalar waterInUpwindCell = 0.0;
        if ( phi > 0  || !mesh_.isInternalFace(fLabel) )
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
        dVf = min(alphaf*phi*dt,waterInUpwindCell);
        isoDebug(Info << "Warning: mag(n0) = " << mag(n0) << " so timeIntegratedFlux calculates dVf from upwind cell alpha value: " << alphaf << endl;)
        return dVf;
    }

    //Find sorted list of times where the isoFace will arrive at face points given initial position x0 and velocity Un0*n0
    const labelList& pLabels = mesh_.faces()[fLabel];
    const label nPoints = pLabels.size();
    pointField fPts(nPoints);
    forAll(fPts,pi)
    {
        fPts[pi] = mesh_.points()[pLabels[pi]];
    }
    scalarField pTimes(fPts.size());
    if ( mag(Un0) > 1e-12 )
    {
        pTimes = ((fPts - x0) & n0)/Un0; //Here we estimate time of arrival to the face points from their normal distance to the initial surface and the surface normal velocity
        dVf = phi/magSf*timeIntegratedArea(fPts,pTimes,dt,magSf,Un0);
       return dVf;
    }
    else //Un0 is almost zero and isoFace is treated as stationary
    {
//        scalar alphaf;
        isoCutFace_.calcSubFace(fLabel,f0);
        scalar alphaf = mag(isoCutFace_.subFaceArea()/magSf);
        dVf = phi*dt*alphaf;
        isoDebug(Info << "Warning: Un0 is almost zero (" << Un0 << ") so calculating dVf on face " << fLabel << " using subFaceFraction giving alphaf = " << alphaf << endl;)
        return dVf;
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
    isoDebug(Info << "Enter timeIntegratedArea for face " << fLabel << endl;)
    scalar tIntArea = 0.0;

    //Sorting face vertex encounter time list
    List<scalar> sortedTimes(pTimes);
    sort(sortedTimes);
    isoDebug(Info << "sortedTimes = " << sortedTimes << endl;)

    //Dealing with case where face is not cut by surface during time interval [0,dt] because face was already passed by surface
    if ( sortedTimes.last() <= 0.0 ) //All cuttings in the past
    {
        isoDebug(Info << "All cuttings in the past" << endl;)
        tIntArea = magSf*dt*pos(Un0); //If all face cuttings were in the past and cell is filling up (Un0>0) then face must be full during whole time interval
        return tIntArea;
    }

    //Dealing with case where face is not cut by surface during time interval [0,dt] because dt is too small for surface to reach closest face point
    if ( sortedTimes.first() >= dt ) //All cuttings in the future
    {
        isoDebug(Info << "All cuttings in the future" << endl;)
        tIntArea = magSf*dt*(1-pos(Un0)); //If all cuttings are in the future but non of them within [0,dt] then if cell is filling up (Un0 > 0) face must be empty during whole time interval
        return tIntArea;
    }

    //Cutting sortedTimes at 0 and dt and sorting out duplicates
    DynamicList<scalar> t(sortedTimes.size()+2);
    t.append(0);
    scalar smallTime = max(1e-6*dt,10*SMALL);
    forAll(sortedTimes,ti)
    {
        if 
        (
            smallTime < sortedTimes[ti] 
                && sortedTimes[ti] < dt-smallTime
                && mag(sortedTimes[ti] - t.last()) > smallTime
        )
        {
            t.append(sortedTimes[ti]);
        }
    }
    t.append(dt);
    isoDebug(Info << "Cutting sortedTimes at 0 and dt: t = " << t << endl;)
//    Info << "times = " << t << endl;
    bool faceUncutInFirstInterval(sortedTimes.first() > 0.0);
    bool faceUncutInLastInterval(sortedTimes.last() < dt);

    //Dealing with cases where face is cut at least once during time interval [0,dt]
    label nt = 0; //sub time interval counter
    scalar initialArea(0.0); //Submerged area at time
    if ( faceUncutInFirstInterval ) //Special treatment for first time interval if face is uncut during this
    {
        tIntArea = magSf*(t[nt+1]-t[nt])*(1.0-pos(Un0)); //If Un0 > 0 cell is filling up - hence if face is cut at a later time but not initially it must be initially empty
        initialArea = magSf*(1.0-pos(Un0));
        isoDebug(Info << "faceUncutInFirstInterval, so special treatment for first time interval: [" << t[nt] << ", " << t[nt+1] << "] giving tIntArea = " << tIntArea << endl;)
        nt++;
    }
    else //calculate initialArea if face is initially cut
    {
        isoDebug(Info << "face is initially cut, so finding initial area, pTimes = " << pTimes << ", Un0 = " << Un0 << endl;)
        isoCutFace_.calcSubFace(fPts,-sign(Un0)*pTimes, 0.0);
        initialArea = mag(isoCutFace_.subFaceArea());
    }
    isoDebug(Info << "InitialArea for next time step corresponds to face phase fraction a0 = " << initialArea/magSf << " where |Sf| = " << magSf << " was used." << endl;)
    while ( nt < t.size()-(1+faceUncutInLastInterval) ) //
    {
        DynamicList<point> cutPoints1(3), cutPoints2(3);
        isoCutFace_.cutPoints(fPts, pTimes, t[nt], cutPoints1);
        isoCutFace_.cutPoints(fPts, pTimes, t[nt+1], cutPoints2);        
        scalar quadArea, intQuadArea;
        quadAreaCoeffs(cutPoints1, cutPoints2, quadArea, intQuadArea);
        scalar integratedQuadArea = sign(Un0)*intQuadArea; //Integration of area(t) = A*t^2+B*t from t = 0 to 1.
        tIntArea += (t[nt+1]-t[nt])*(initialArea + integratedQuadArea);
        initialArea += sign(Un0)*quadArea; //Adding quad area to submerged area
        isoDebug(Info << "Integrating area for " << nt+1 << "'th time interval: [" << t[nt] << ", " << t[nt+1] << "] giving tIntArea = " << tIntArea << " and a0 = " << initialArea/magSf << endl;)
        nt++;
    }

    if ( faceUncutInLastInterval && (nt+1) < t.size() ) //Special treatment for last time interval if face is uncut during this
    {
        isoDebug(Info << "faceUncutInLastInterval, so special treatment for last (" << nt+1 << "'th) time interval: [" << t[nt] << ", " << t[nt+1] << "]" << endl;)
        tIntArea += magSf*(t[nt+1]-t[nt])*pos(Un0); //If face is cut at some intermediate time but not at last time, then if Un0 > 0 (cell filling up) face must be filled at last time interval.
    }
    return tIntArea;
}


void Foam::isoAdvection::getDownwindFaces
(
    const label ci,
    DynamicList<label>& downwindFaces
)
{
    isoDebug(Info << "Enter getDownwindFaces " << endl;)
    const cellList& cells = mesh_.cells();
    const labelList& fLabels = cells[ci];

    // Check all faces of the cell
    forAll(fLabels,fi)
    {
        const label fLabel = fLabels[fi];
        const scalar phi = faceValue(phi_,fLabel);

        if (mesh_.owner()[fLabel] == ci)
        {
//            if (phi > 1e-12)
            if (phi > 10*SMALL)
            {
                downwindFaces.append(fLabel);
            }
        }
//        else if ( phi < -1e-12 ) //ci must be neighbour of fLabel
        else if ( phi < -10*SMALL ) //ci must be neighbour of fLabel
        {
            downwindFaces.append(fLabel);
        }
    }
    downwindFaces.shrink();
}


bool Foam::isoAdvection::isADownwindFace
(
    const label faceI,
    const label cellI
)
{
    const scalar phi = faceValue(phi_,faceI);

    bool isDownwindFace = false;
    if (mesh_.owner()[faceI] == cellI)
    {
//        if ( phi > 1e-12 )
        if ( phi > 10*SMALL )
        {
            isDownwindFace = true;
        }
    }
//    else if ( phi < -1e-12 ) //cellI assumed to be neighbour of faceI
    else if ( phi < -10*SMALL ) //cellI assumed to be neighbour of faceI
    {
        isDownwindFace = true;
    }

    return isDownwindFace;
}


void Foam::isoAdvection::quadAreaCoeffs
(
    const DynamicList<point>& pf0,
    const DynamicList<point>& pf1,
    scalar& quadArea,
    scalar& intQuadArea
)
{
    isoDebug(Info << "Enter quadAreaCoeffs" << endl;)

    label np0(pf0.size()), np1(pf1.size());
//    isoDebug(Info << "Face " << fLabel << " was cut at " << pf0 << " by f0 = " << f0 << " and at " << pf1 << " by " << " f1 = " << f1 << endl;)

    scalar alpha = 0.0;
    scalar beta = 0.0;
    quadArea = 0.0;
    intQuadArea = 0.0;

    if ( np0 > 0 && np1 > 0)
    {
        //Defining quadrilateral vertices
        vector A(pf0[0]), C(pf1[0]), B(vector::zero), D(vector::zero);

        //Triangle cases
        if (np0 == 2 && mag(pf0[0]-pf0[1]) > SMALL)
        {
            B = pf0[1];
        }
        else
        {
            B = A + 1e-4*(pf1[1]-pf1[0]);
            if ( np0  != 1 )
            {
                Info << "Warning: Vertex face was cut at pf0 = " << pf0 << endl;
            }
        }

        if (np1 == 2 && mag(pf1[0]-pf1[1]) > SMALL)
        {
            D = pf1[1];
        }
        else
        {
            D = C + 1e-4*(A-B);
            if ( np1 != 1)
            {
                Info << "Warning: Vertex face was cut at pf1 = " << pf1 << endl;
            }
        }
        
        //Defining local coordinates for area integral calculation
        vector xhat = B-A;
        xhat /= mag(xhat);
        vector yhat = (D-A);
        yhat -= ((yhat & xhat)*xhat);
        yhat /= mag(yhat);
        
//        isoDebug(Info << "xhat = " << xhat << ", yhat = " << yhat << ", zhat = " << zhat << ". x.x = " << (xhat & xhat) << ", y.y = " << (yhat & yhat) <<", z.z = " << (zhat & zhat) << ", x.y = " << (xhat & yhat) << ", x.z = " << (xhat & zhat) << ", y.z = " << (yhat & zhat) << endl;)


        //Swapping pf1 points if pf0 and pf1 point in same general direction (because we want a quadrilateral ABCD where pf0 = AB and pf1 = CD)
        if ( ((B-A) & (D-C)) > 0 )
        {
////            Info << "Swapping C and D" << endl;
            vector tmp = D;
            D = C;
            C = tmp;
        }

////        Info << "A = " << A << ", B = " << B << ", C = " << C << ", D = " << D << endl;
        scalar Bx = mag(B-A);
        scalar Cx = (C-A) & xhat;
        scalar Cy = mag((C-A) & yhat);
        scalar Dx = (D-A) & xhat;
        scalar Dy = mag((D-A) & yhat);


//      area = ((Cx-Bx)*Dy-Dx*Cy)/6.0 + 0.25*Bx*(Dy+Cy);
        alpha = 0.5*((Cx-Bx)*Dy-Dx*Cy);
        beta = 0.5*Bx*(Dy+Cy);
        quadArea = alpha + beta;
        intQuadArea = (1/3)*alpha + 0.5*beta;

////        Info << "Bx = " << Bx << ", Cx = " << Cx << ", Cy = " << Cy << ", Dx = " << Dx << ", Dy = " << Dy << ", alpha = " << alpha << ", beta = " << beta << endl;
        //area(t) = A*t^2+B*t
        //integratedArea = A/3+B/2
    }
    else
    {
        Info << "Vertex face was cut at " << pf0
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


void Foam::isoAdvection::limitFluxes
(
    surfaceScalarField& dVf,
    const scalar dt
)
{
//    scalarField alphaNew = alpha1In_ - fvc::surfaceIntegrate(dVf);
    scalar aTol = 1.0e-12;
    scalar maxAlphaMinus1 = 1;//max(alphaNew-1);
    scalar minAlpha = -1;//min(alphaNew);
    label nUndershoots = 20;//sum(neg(alphaNew+aTol));
    label nOvershoots = 20;//sum(pos(alphaNew-1-aTol));
    cellIsBounded_ = false;
    
    for (label n = 0; n < nAlphaBounds_; n++)
    {
        Info << "Running bounding number " << n + 1 << " of time " 
            << mesh_.time().value() << endl;

        if ( maxAlphaMinus1 > aTol )
        {
            isoDebug(Info << "Bound from above... " << endl;)
//          scalarField dVfcorrected = dVf.internalField();
            surfaceScalarField dVfcorrected("dVfcorrected", dVf);
            DynamicList<label> correctedFaces(3*nOvershoots);
            boundFromAbove(alpha1In_,dt,dVfcorrected,correctedFaces);
            forAll(correctedFaces,fi)
            {
                label fLabel = correctedFaces[fi];
                //Change to treat boundaries consistently
                faceValue(dVf,fLabel,faceValue(dVfcorrected,fLabel));
            }
            syncProcPatches(dVf,phi_);
        }

        if ( minAlpha < -aTol )
        {
            isoDebug(Info << "Bound from below... " << endl;)
            scalarField alpha2 = 1.0 - alpha1In_;
            surfaceScalarField dVfcorrected("dVfcorrected", phi_*dimensionedScalar("dt", dimTime, dt) - dVf);
//          dVfcorrected -= dVf; //phi_ and dVf have same sign and dVf is the portion of phi_*dt that is water.
            //If phi_ > 0 then dVf > 0 and mag(phi_*dt-dVf) < mag(phi_*dt) as it should.
            //If phi_ < 0 then dVf < 0 and mag(phi_*dt-dVf) < mag(phi_*dt) as it should.
            DynamicList<label> correctedFaces(3*nUndershoots);
            boundFromAbove(alpha2,dt,dVfcorrected,correctedFaces);
            forAll(correctedFaces,fi)
            {
                label fLabel = correctedFaces[fi];
                //Change to treat boundaries consistently
                scalar phi = faceValue(phi_,fLabel);
                scalar dVcorr = faceValue(dVfcorrected,fLabel);
                faceValue(dVf,fLabel,phi*dt - dVcorr);
            }
            syncProcPatches(dVf,phi_);
        }
/*
        //Check if still unbounded
        alphaNew = alpha1In_ - fvc::surfaceIntegrate(dVf);
        maxAlphaMinus1 = max(alphaNew-1);
        minAlpha = min(alphaNew);
        nUndershoots = sum(neg(alphaNew+aTol));
        nOvershoots = sum(pos(alphaNew-1-aTol));
        Info << "After bounding number " << n+1 << " of time " << mesh_.time().value() << ":" << endl;
        Info << "nOvershoots = " << nOvershoots << " with max(alphaNew-1) = " << maxAlphaMinus1 << " and nUndershoots = " << nUndershoots << " with min(alphaNew) = " << minAlpha << endl;
*/
    }
}


void Foam::isoAdvection::boundFromAbove
(
    const scalarField& alpha1,
    const scalar dt,
    surfaceScalarField& dVf,
    DynamicList<label>& correctedFaces
)
{
    correctedFaces.clear();
//    scalar aTol = 1e-12;
    scalar aTol = 10*SMALL;

    forAll(alpha1,ci)
    {
        if (checkBounding_[ci])
        {

            scalar Vi = mesh_.V()[ci];
            scalar alpha1New = alpha1[ci] - netFlux(dVf,ci)/Vi;
            scalar alphaOvershoot = alpha1New - 1.0;
            scalar fluidToPassOn = alphaOvershoot*Vi;
            label nFacesToPassFluidThrough = 1;

            //First try to pass surplus fluid on to neighbour cells that are not filled and to which dVf < phi*dt
            while ( alphaOvershoot > aTol && nFacesToPassFluidThrough > 0 )
            {
                cellIsBounded_[ci] = true;
                isoDebug(Info << "\n\nBounding cell " << ci << " with alpha overshooting " << alphaOvershoot << endl;)
                //First find potential neighbour cells to pass surplus water to
                DynamicList<label> downwindFaces(mesh_.cells()[ci].size());
                getDownwindFaces(ci,downwindFaces);
                DynamicList<label> facesToPassFluidThrough(downwindFaces.size());
                DynamicList<scalar> dVfmax(downwindFaces.size());
                DynamicList<scalar> phi(downwindFaces.size());
                scalar dVftot = 0.0;
                nFacesToPassFluidThrough = 0;

                isoDebug(Info << "downwindFaces: " << downwindFaces << endl;)
                forAll(downwindFaces,fi)
                {
                    label fLabel = downwindFaces[fi];
                    scalar maxExtraFaceFluidTrans = 0.0;
                    scalar phif = faceValue(phi_,fLabel);
                    scalar dVff = faceValue(dVf,fLabel);
                    maxExtraFaceFluidTrans = mag(phif*dt - dVff);

                    //dVf has same sign as phi and so if phi>0 we have mag(phi_[fLabel]*dt) - mag(dVf[fLabel]) = phi_[fLabel]*dt - dVf[fLabel]
                    //If phi<0 we have mag(phi_[fLabel]*dt) - mag(dVf[fLabel]) = -phi_[fLabel]*dt - (-dVf[fLabel]) > 0 since mag(dVf) < phi*dt
                    isoDebug(Info << "downwindFace " << fLabel << " has maxExtraFaceFluidTrans = " << maxExtraFaceFluidTrans << endl;)
                    if ( maxExtraFaceFluidTrans/Vi > aTol )
    //              if ( maxExtraFaceFluidTrans/Vi > aTol && mag(dVfIn[fLabel])/Vi > aTol ) //Last condition may be important because without this we will flux through uncut downwind faces
                    {
                        facesToPassFluidThrough.append(fLabel);
                        phi.append(phif);
                        dVfmax.append(maxExtraFaceFluidTrans);
                        dVftot += mag(phif*dt);
                    }
                }

                isoDebug(Info << "\nfacesToPassFluidThrough: " << facesToPassFluidThrough << ", dVftot = " << dVftot << " m3 corresponding to dalpha = " << dVftot/Vi << endl;)
                forAll(facesToPassFluidThrough,fi)
                {
                    const label fLabel = facesToPassFluidThrough[fi];
                    scalar fluidToPassThroughFace = fluidToPassOn*mag(phi[fi]*dt)/dVftot;
                    nFacesToPassFluidThrough += pos(dVfmax[fi] - fluidToPassThroughFace);
                    fluidToPassThroughFace = min(fluidToPassThroughFace,dVfmax[fi]);
                    scalar dVff = faceValue(dVf,fLabel);
                    dVff += sign(phi[fi])*fluidToPassThroughFace;
                    faceValue(dVf,fLabel,dVff);
                    checkIfOnProcPatch(fLabel);
                    correctedFaces.append(fLabel);
                }
                alpha1New = alpha1[ci] - netFlux(dVf,ci)/Vi;
                alphaOvershoot = alpha1New - 1.0;
                fluidToPassOn = alphaOvershoot*Vi;
                isoDebug(Info << "\nNew alpha for cell " << ci << ": " << alpha1New << endl;)
            }
        }
    }
    isoDebug(Info << "correctedFaces = " << correctedFaces << endl;)
}


Foam::scalar Foam::isoAdvection::netFlux
(
    const surfaceScalarField& dVf,
    const label cLabel
)
{
    scalar dV = 0.0;

    const labelList& fLabels = mesh_.cells()[cLabel];

    forAll (fLabels, fi)
    {
        const label fLabel = fLabels[fi];
        const scalar dVff = faceValue(dVf,fLabel);

        if (mesh_.owner()[fLabel] == cLabel)
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


Foam::scalar Foam::isoAdvection::faceValue
(
    const surfaceScalarField& f,
    const label fLabel
)
{
    isoDebug(Info << "Enter faceValue reading scalar " << endl;)
    if (mesh_.isInternalFace(fLabel))
    {
        return scalar(f.internalField()[fLabel]);
    }
    else
    {
        // Boundary face.  Find out which face of which patch
        const label patchI = mesh_.boundaryMesh().whichPatch(fLabel);

        // Handle empty patches
        if (mesh_.boundary()[patchI].size() == 0)
        {
            return 0;
        }

        if (patchI < 0 || patchI >= mesh_.boundaryMesh().size())
        {
            FatalErrorIn
            (
                "void isoAdvection::faceValue(...)"
            )   << "Cannot find patch for face " << fLabel
                << abort(FatalError);
        }

        const label faceI =
            mesh_.boundaryMesh()[patchI].whichFace
            (
                fLabel
            );

        return f.boundaryField()[patchI][faceI];
    }
}


void Foam::isoAdvection::faceValue
(
    surfaceScalarField& f,
    const label fLabel,
    const scalar value
)
{
    isoDebug(Info << "Enter faceValue setting scalar" << endl;)
    if (mesh_.isInternalFace(fLabel))
    {
        f.internalField()[fLabel] = value;
    }
    else
    {
        // Boundary face.  Find out which face of which patch
        const label patchI = mesh_.boundaryMesh().whichPatch(fLabel);

        // Handle empty patches
        if (mesh_.boundary()[patchI].size() == 0)
        {
            return;
        }

        if (patchI < 0 || patchI >= mesh_.boundaryMesh().size())
        {
            FatalErrorIn
            (
                "void isoAdvection::faceValue(...)"
            )   << "Cannot find patch for face " << fLabel
                << abort(FatalError);
        }

        const label faceI =
            mesh_.boundaryMesh()[patchI].whichFace
            (
                fLabel
            );

        f.boundaryField()[patchI][faceI] = value;
    }
}


Foam::vector Foam::isoAdvection::faceValue
(
    const surfaceVectorField& f,
    const label fLabel
)
{
    isoDebug(Info << "Enter faceValue reading vector" << endl;)
    if (mesh_.isInternalFace(fLabel))
    {
        return vector(f.internalField()[fLabel]);
    }
    else
    {
        // Boundary face.  Find out which face of which patch
        const label patchI = mesh_.boundaryMesh().whichPatch(fLabel);

        // Handle empty patches
        if (mesh_.boundary()[patchI].size() == 0)
        {
            return vector::zero;
        }

        if (patchI < 0 || patchI >= mesh_.boundaryMesh().size())
        {
            FatalErrorIn
            (
                "void isoAdvection::faceValue(...)"
            )   << "Cannot find patch for face " << fLabel
                << abort(FatalError);
        }

        const label faceI =
            mesh_.boundaryMesh()[patchI].whichFace
            (
                fLabel
            );

        return f.boundaryField()[patchI][faceI];
    }
}


void Foam::isoAdvection::faceValue
(
    surfaceVectorField& f,
    const label fLabel,
    const vector value
)
{
    isoDebug(Info << "Enter faceValue setting vector" << endl;)
    if (mesh_.isInternalFace(fLabel))
    {
        f.internalField()[fLabel] = value;
    }
    else
    {
        // Boundary face.  Find out which face of which patch
        const label patchI = mesh_.boundaryMesh().whichPatch(fLabel);

        // Handle empty patches
        if (mesh_.boundary()[patchI].size() == 0)
        {
            return;
        }

        if (patchI < 0 || patchI >= mesh_.boundaryMesh().size())
        {
            FatalErrorIn
            (
                "void isoAdvection::faceValue(...)"
            )   << "Cannot find patch for face " << fLabel
                << abort(FatalError);
        }

        const label faceI =
            mesh_.boundaryMesh()[patchI].whichFace
            (
                fLabel
            );

        f.boundaryField()[patchI][faceI] = value;
    }
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
        forAll(procPatchLabels_, patchLabelI)
        {
            const label patchI = procPatchLabels_[patchLabelI];

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchI]);

            UOPstream toNbr(procPatch.neighbProcNo(), pBufs);
            const scalarField& pFlux = dVf.boundaryField()[patchI];

            const List<label>& surfCellFacesOnProcPatch = surfaceCellFacesOnProcPatches_[patchI];
            List<scalar> dVfPatch(surfCellFacesOnProcPatch.size());
            forAll(dVfPatch,faceI)
            {
                dVfPatch[faceI] = pFlux[surfCellFacesOnProcPatch[faceI]];
            }

//          Pout << "Sent at time = " << mesh_.time().value() << ": surfCellFacesOnProcPatch = " << surfCellFacesOnProcPatch << endl;
//          Pout << "Sent at time = " << mesh_.time().value() << ": dVfPatch = " << dVfPatch << endl;

/*
            toNbr << pFlux
                  << surfCellFacesOnProcPatch
                  << dVfPatch;
*/
            toNbr << surfCellFacesOnProcPatch
                  << dVfPatch;
        }

        pBufs.finishedSends();


        // Receive and combine.
        forAll(procPatchLabels_, patchLabelI)
        {
            const label patchI = procPatchLabels_[patchLabelI];

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchI]);

            UIPstream fromNeighb(procPatch.neighbProcNo(), pBufs);
            DynamicList<label> fLabels(100);
            DynamicList<scalar> nbrdVfs(100);
//          const scalarField nbrFlux(fromNeighb);

/*
            scalarField nbrFlux(procPatch.size());

            fromNeighb >> nbrFlux
                       >> fLabels
                       >> nbrdVfs;
*/

            fromNeighb >> fLabels
                       >> nbrdVfs;

//          Pout << "Recieved at time = " << mesh_.time().value() << ": surfCellFacesOnProcPatch = " << fLabels << endl;
//          Pout << "Recieved at time = " << mesh_.time().value() << ": dVfPatch = " << nbrdVfs << endl;

            // Combine fluxes
            scalarField& localFlux = dVf.boundaryField()[patchI];
/*

            const scalarField& phib = phi.boundaryField()[patchI];
            forAll (nbrFlux, faceI)
            {
                if (phib[faceI] < 0)
                {
                    localFlux[faceI] = -nbrFlux[faceI];
                }
            }
*/
            forAll (fLabels, faceI)
            {
                const label fLabel = fLabels[faceI];
                localFlux[fLabel] = - nbrdVfs[faceI];
                if (mag(localFlux[fLabel] + nbrdVfs[faceI]) > 10*SMALL)
                {
                    Pout << "localFlux[fLabel] = " << localFlux[fLabel] << " and nbrdVfs[faceI] = " << nbrdVfs[faceI] << " for fLabel = " << fLabel << endl;
                }
            }


        }
/*
        //Write out results for checking
        forAll(procPatchLabels_, patchLabelI)
        {
            const label patchI = procPatchLabels_[patchLabelI];
            scalarField& localFlux = dVf.boundaryField()[patchI];
            Pout << "time = " << mesh_.time().value() << ": localFlux = " << localFlux << endl;

        }
*/
        //Reinitialising list used for minimal parallel communication
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
   if ( !mesh_.isInternalFace(faceI) )
    {
        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        const label patchI = patches.whichPatch(faceI);

        if
        (
            isA<processorPolyPatch>(patches[patchI])
         && patches[patchI].size() > 0
        )
        {
            const label fLabel =
            patches[patchI].whichFace
            (
                faceI
            );
//          Pout << "fLabel = " << fLabel << " for faceI = " << faceI << " on patchI = " << patchI << endl;
            surfaceCellFacesOnProcPatches_[patchI].append(fLabel);
        }
    }
}


void Foam::isoAdvection::getTransportedVolume
(
    const scalar dt,
    surfaceScalarField& dVf
)
{
    isoDebug(Info << "Enter getTransportedVolume" << endl;)

    //Interpolating alpha1 cell centre values to mesh points (vertices)
    ap_ = vpi_.interpolate(alpha1_);

    //Initialising dVf with upwind values, i.e. phi[fLabel]*alpha1[upwindCell]*dt    
    dVf = upwind<scalar>(mesh_, phi_).flux(alpha1_)*dimensionedScalar("dt", dimTime, dt);

    //Do the isoAdvection on surface cells
    timeIntegratedFlux(dt, dVf);

    //Syncronize processor patches
    syncProcPatches(dVf,phi_);

    //Adjust dVf for unbounded cells
    limitFluxes(dVf,dt);
}


void Foam::isoAdvection::advect
(
    const scalar dt
)
{
    isoDebug(Info << "Enter advect" << endl;)

    //Interpolating alpha1 cell centre values to mesh points (vertices)
    ap_ = vpi_.interpolate(alpha1_);

    //Initialising dVf with upwind values, i.e. phi[fLabel]*alpha1[upwindCell]*dt    
    dVf_ = upwind<scalar>(mesh_, phi_).flux(alpha1_)*dimensionedScalar("dt", dimTime, dt);

    //Do the isoAdvection on surface cells
    timeIntegratedFlux(dt, dVf_);

    //Syncronize processor patches
    syncProcPatches(dVf_,phi_);

    //Adjust dVf for unbounded cells
    limitFluxes(dVf_,dt);
    
    alpha1_ -= fvc::surfaceIntegrate(dVf_);
    alpha1_.correctBoundaryConditions();
}


void Foam::isoAdvection::getSurfaceCells
(
    cellSet& surfCells
)
{
    surfCells.clear();
    forAll(surfCells_,i)
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
    forAll(cellIsBounded_,i)
    {
        if (cellIsBounded_[i])
        {
            boundCells.insert(i);
        }
    }
}


void Foam::isoAdvection::getRhoPhi
(
    surfaceScalarField& rhoPhi,
    const dimensionedScalar rho1,
    const dimensionedScalar rho2,
    const dimensionedScalar dt
)
{
    rhoPhi = (rho1 - rho2)*dVf_/dt + rho2*phi_;
}

// ************************************************************************* //
