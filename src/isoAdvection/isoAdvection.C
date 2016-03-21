/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "isoAdvection.H"
#include "isoCutter.H"
#include "volPointInterpolation.H"
#include "interpolationCellPoint.H"
#include "fvcSurfaceIntegrate.H"

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
    const volVectorField& U,
    const dictionary& dict
)
:
    mesh_(alpha1.mesh()),
    alpha1_(alpha1),
    vpi_(mesh_),
    ap_(mesh_.nPoints(),0.0),
    phi_(phi),
    U_(U),
    isSurfaceCell_(mesh_.nCells(),false),
    nAlphaBounds_(dict.lookupOrDefault<label>("nAlphaBounds", 1)),
    vof2IsoTol_(dict.lookupOrDefault<scalar>("vof2IsoTol", 1e-8)),
    surfCellTol_(dict.lookupOrDefault<scalar>("surfCellTol", 1e-8)),
    writeToLog_(dict.lookupOrDefault<bool>("writeToLog", true))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::isoAdvection::timeIntegratedFlux
(
    const scalar dt,
    surfaceScalarField& dVf
)
{
    isoDebug(Info << "Enter timeIntegratedFlux" << endl;)

    //Interpolating VOF field to mesh points
/*
    //Testing
    const labelListList& pCells = mesh_.pointCells();
    forAll(pCells,pi)
    {
        const labelList& pci = pCells[pi];
        ap_[pi] = 0.0;
        scalar vol(0.0);
        forAll(pci,ci)
        {
            const label cLabel = pci[ci];
            ap_[pi] += alpha1_[cLabel]*mesh_.V()[cLabel];
            vol += mesh_.V()[cLabel];
        }
        ap_[pi] /= vol;
    }
//  isoDebug(Info << "Using cell volume weighted cell-point interpolation" << endl;)
*/

//    volPointInterpolation vpi(mesh_);
    ap_ = vpi_.interpolate(alpha1_);

    //Construction of isosurface calculator to get access to its functionality
    isoCutter cutter(mesh_,ap_);

    //Find surface cells
    DynamicList<label> surfaceCells(ceil(mesh_.nCells()/10));
    findSurfaceCells(surfaceCells);

    //Estimated total water volume transported across mesh faces during time interval dt. The sign convenctino is like the flux phi, i.e. positive means out of owner cell.
    scalarField& dVfIn = dVf.internalField();
    const scalarField& phiIn = phi_.internalField();
    const scalarField& alphaIn = alpha1_.internalField();

    const unallocLabelList& own = mesh_.owner();
    const unallocLabelList& nei = mesh_.neighbour();

    //Upwinding on all internal faces receiving fluid from a non-surface cell
    forAll (dVfIn, faceI)
    {
        dVfIn[faceI] = 0;
        label donorCell(-1);

        if ( phiIn[faceI] >= 1e-12 )
        {
            donorCell = own[faceI];
        }
        else
        {
            donorCell = nei[faceI];
        }

        dVf[faceI] = phiIn[faceI]*max(min(1, alphaIn[donorCell]), 0)*dt;
    }

    forAll (dVf.boundaryField(), patchI)
    {
        if (mesh_.boundary()[patchI].size() > 0)
        {
            fvsPatchScalarField& patchdVF = dVf.boundaryField()[patchI];
            const fvsPatchScalarField& patchPhi = phi_.boundaryField()[patchI];
            const fvPatchScalarField& patchAlpha1 = alpha1_.boundaryField()[patchI];

            const unallocLabelList& fc = mesh_.boundary()[patchI].faceCells();

            forAll (patchdVF, pFaceI)
            {
                if (patchPhi[pFaceI] >= 0)
                {
                    patchdVF[pFaceI] = patchPhi[pFaceI]*
                        Foam::max(Foam::min(scalar(1), alphaIn[fc[pFaceI]]), scalar(0))*dt;
                }
                else
                {
                    patchdVF[pFaceI] = patchPhi[pFaceI]*
                        Foam::max(Foam::min(scalar(1), patchAlpha1[pFaceI]), scalar(0))*dt;
                }
            }
        }
    }

    //Isocutting estimate of water flux from surface cellss
    interpolationCellPoint<vector> UInterp(U_);

    forAll(surfaceCells,cellI)
    {
        const label ci = surfaceCells[cellI];
        isoDebug(Info << "\n------------ Cell " << ci << " with alpha1 = " << alpha1_[ci] << " and 1-alpha1 = " << 1.0-alpha1_[ci] << " ------------" << endl;)

        //Make list of all cell faces out of which fluid is flowing
        DynamicList<label> downwindFaces((mesh_.cells()[ci]).size());
        getDownwindFaces(ci, downwindFaces);
        isoDebug(Info << "downwindFaces: " << downwindFaces << "\n" << endl;)

        //Calculate isoFace0 centre xs0, normal ns0, and velocity U0 = U(xs0)
        scalar f0(0.0), Un0(0.0);
        vector x0(vector::zero), n0(vector::zero);
        calcIsoFace(ci,x0,n0,f0,Un0,UInterp); //This one really also should give us a0 on all faces since it is calculated anyway. Do this with a cutCell structure
        isoDebug(Info << "calcIsoFace gives initial surface: \nx0 = " << x0 << ", \nn0 = " << n0 << ", \nf0 = " << f0 << ", \nUn0 = " << Un0 << endl;)

        //Estimating time integrated water flux through each downwinding face
        forAll(downwindFaces,fi)
        {
            const label fLabel = downwindFaces[fi];
            isoDebug(Info << "\nFace " << fLabel << " with outward normal nf = " << sign( (mesh_.Cf()[fLabel]-mesh_.C()[ci]) & mesh_.Sf()[fLabel] )*mesh_.Sf()[fLabel]/mag(mesh_.Sf()[fLabel]) << endl;)
            const scalar tif = timeIntegratedFlux(fLabel,x0,n0,Un0,f0,dt);
            faceValue(dVf,fLabel,tif);
        }
    }
}


void Foam::isoAdvection::findSurfaceCells
(
    DynamicList<label>& surfaceCells
)
{
    isoDebug(Info << "Enter findSurfaceCells" << endl;)

    isSurfaceCell_ = false;
    forAll(alpha1_,ci)
    {
        if ( ( surfCellTol_ < alpha1_[ci] && alpha1_[ci] < 1.0-surfCellTol_ ) )
        {
            isSurfaceCell_[ci] = true;
            surfaceCells.append(ci);
        }
    }
    surfaceCells.shrink();
    Info << "\nnSurfaceCells = " << surfaceCells.size() << endl;
}


void Foam::isoAdvection::calcIsoFace
(
    const label ci,
    vector& x0,
    vector& n0,
    scalar& f0,
    scalar& Un0,
    const interpolationCellPoint<vector>& UInterp
)
{
    isoDebug(Info << "Enter calcIsoFace" << endl;)
    //Construction of isosurface calculator to get access to its functionality
    isoCutter cutter(mesh_,ap_);
    label maxIter(100);
    vector subCellCtr;
    cutter.vofCutCell(ci, alpha1_[ci], vof2IsoTol_, maxIter, f0, subCellCtr);
    cutter.isoFaceCentreAndArea(ci,f0,x0,n0); //Stupid to recalculate this here - should be provided by vofCutCell above
    isoDebug(Info << "f0 = " << f0 << ", x0 = " << x0 << ", n0 = " << n0 << endl;)

    if ( mag(n0) < 1.0e-8 ) //Cell almost full or empty so isoFace is undefined. Calculating normal by going a little into the cell
    {
        scalar aMin(GREAT), aMax(-GREAT);
        const labelList cellPts = mesh_.cellPoints()[ci];
        subSetExtrema(ap_,cellPts,aMin,aMax);
        if (mag(f0-aMin) < mag(f0-aMax)) //f0 almost equal to aMin, i.e. cell almost full
        {
            isoDebug(Info << "Cell is almost full with aMin = " << aMin << " and aMax = " << aMax << endl;)
            scalar f0Inside =  aMin + 1e-3*(aMax-aMin);
            vector xdum;
            cutter.isoFaceCentreAndArea(ci,f0Inside,xdum,n0);
            if (((x0 - mesh_.C()[ci]) & n0) < 0.0) //n0 should point from cell centre towards isoFace centre for an almost full cell
            {
                n0 *= (-1.0);
                isoDebug(Info << "Changing direction of n0 " << n0 << endl;)
            }
        }
        else if (mag(f0-aMin) > mag(f0-aMax))//f0 almost equal to aMax, i.e. cell almost empty
        {
            isoDebug(Info << "Cell is almost empty with aMin = " << aMin << " and aMax = " << aMax << endl;)
            scalar f0Inside =  aMax + 1e-3*(aMin-aMax);
            vector xdum;
            cutter.isoFaceCentreAndArea(ci,f0Inside,xdum,n0);
            if (((mesh_.C()[ci] - x0) & n0) < 0.0) //n0 should point from isoFace centre towards cell centre for an almost empty cell
            {
                n0 *= (-1.0);
                isoDebug(Info << "Changing direction of n0 " << n0 << endl;)
            }
        }
        else
        {
            Info << "Warning: aMax = aMin = " << aMax << " for cell " << ci << endl;
        }

    } //Make sure n0 points out of subCell
    else if (((x0 - subCellCtr) & n0) < 0.0) //n0 should point out of water surface, i.e. in the direction from subCell centre to isoFace centre
    {
        n0 *= (-1.0);
        isoDebug(Info << "Changing direction of n0 " << n0 << endl;)
    }

    if ( mag(n0) > 1.0e-8 )
    {
        isoDebug(Info << "Normalising n0: " << n0 << endl;)
        n0 /= mag(n0);
    }

    //Interpolate velocity to isoFace centre
    vector U0 = UInterp.interpolate(x0,ci);
    Un0 = (U0 & n0);
}


Foam::scalar Foam::isoAdvection::timeIntegratedFlux
(
    const label fLabel,
    const vector& x0,
    const vector& n0,
    const scalar Un0,
    const scalar f0,
    const scalar dt
)
{
    isoDebug(Info << "Enter timeIntegratedFlux for face " << fLabel << endl;)

    scalar dVf(0.0); //Volume flowing through face in time interval [0,dt] to be calculated below
    const scalar phi = faceValue(phi_,fLabel);

    if (mag(n0) < .5)
    {
        scalar alphaf = 0.0;
        scalar waterInUpwindCell = 0.0;
        if ( phi > 0  || !mesh_.isInternalFace(fLabel) )
        {
            label upwindCell = mesh_.owner()[fLabel];
            alphaf = alpha1_[upwindCell];
            waterInUpwindCell = alphaf*mesh_.V()[upwindCell];
        }
        else
        {
            label upwindCell = mesh_.neighbour()[fLabel];
            alphaf = alpha1_[upwindCell];
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
        dVf = phi/mag(mesh_.Sf()[fLabel])*timeIntegratedArea(fLabel,fPts,pTimes,dt,Un0);
        return dVf;
    }
    else //Un0 is almost zero
    {
        isoCutter cutter(mesh_,ap_);
        scalar alphaf;
        cutter.subFaceFraction(fLabel,f0,alphaf);
        dVf = phi*dt*alphaf;
        isoDebug(Info << "Warning: Un0 is almost zero (" << Un0 << ") so calculating dVf on face " << fLabel << " using subFaceFraction giving alphaf = " << alphaf << endl;)
        return dVf;
    }
}


Foam::scalar Foam::isoAdvection::timeIntegratedArea
(
    const label fLabel,
    const pointField& fPts,
    const scalarField& pTimes,
    const scalar dt,
    const scalar Un0
)
{
    isoDebug(Info << "Enter timeIntegratedArea for face " << fLabel << endl;)
    scalar tIntArea = 0.0;
    const label nPoints = fPts.size();

    //Face area
    scalar magSf;
    if ( (mesh_.faces()[fLabel]).size() == fPts.size() ) //full face
    {
        magSf = mag(mesh_.Sf()[fLabel]);
    }
    else //triangular subface
    {
        magSf = mag((fPts[1]-fPts[0])^(fPts[2]-fPts[0]));
    }

    //Sorting face vertex encounter time list
    scalarField sortedTimes(pTimes);
    sort(sortedTimes);
    isoDebug(Info << "sortedTimes = " << sortedTimes << endl;)

    //Dealing with case where face is not cut by surface during time interval [0,dt] because face was already passed by surface
    if ( sortedTimes[nPoints-1] <= 0.0 ) //All cuttings in the past
    {
        isoDebug(Info << "All cuttings in the past" << endl;)
        tIntArea = magSf*dt*pos(Un0); //If all face cuttings were in the past and cell is filling up (Un0>0) then face must be full during whole time interval
        return tIntArea;
    }

    //Dealing with case where face is not cut by surface during time interval [0,dt] because dt is too small for surface to reach closest face point
    if ( sortedTimes[0] >= dt ) //All cuttings in the future
    {
        isoDebug(Info << "All cuttings in the future" << endl;)
        tIntArea = magSf*dt*(1-pos(Un0)); //If all cuttings are in the future but non of them within [0,dt] then if cell is filling up (Un0 > 0) face must be empty during whole time interval
        return tIntArea;
    }

    //Cutting sortedTimes at 0 and dt
    DynamicList<scalar> t(sortedTimes.size()+2);
    t.append(0.0);
    forAll(sortedTimes,ti)
    {
        if ( 1e-3*dt < sortedTimes[ti] && sortedTimes[ti] < (1-1e-3)*dt )
        {
            t.append(sortedTimes[ti]);
        }
    }
    t.append(dt);
    isoDebug(Info << "Cutting sortedTimes at 0 and dt: t = " << t << endl;)

    bool faceUncutInFirstInterval(sortedTimes[0] > 0.0);
    bool faceUncutInLastInterval(sortedTimes[nPoints-1] < dt);

    //Dealing with cases where face is cut at least once during time interval [0,dt]
    label nt = 0; //Time interval counter
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
        isoCutter cutter(mesh_,ap_);
        initialArea = cutter.getSubFaceFraction(fLabel, -sign(Un0)*pTimes, 0.0)*magSf; //does not use ap_
    }

    isoDebug(Info << "InitialArea for next time step corresponds to face phase fraction a0 = " << initialArea/magSf << " where |Sf| = " << magSf << " was used." << endl;)
    while ( nt < t.size()-(1+faceUncutInLastInterval) ) //
    {
        scalar quadArea(0.0), intQuadArea(0.0);
        getQuadArea(fLabel,fPts,pTimes,t[nt],t[nt+1],quadArea,intQuadArea);
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
    const cellList& cells = mesh_.cells();
    const labelList& fLabels = cells[ci];

    // Check all faces of the cell
    forAll(fLabels,fi)
    {
        const label fLabel = fLabels[fi];
        const scalar phi = faceValue(phi_,fLabel);

        if (mesh_.owner()[fLabel] == ci)
        {
            if (phi > 1e-12)
            {
                downwindFaces.append(fLabel);
            }
        }
        else if ( phi < -1e-12 ) //ci must be neighbour of fLabel
        {
            downwindFaces.append(fLabel);
        }
    }
    downwindFaces.shrink();
}


void Foam::isoAdvection::getQuadArea
(
    const label fLabel,
    const pointField& p,
    const scalarField& f,
    const scalar f0,
    const scalar f1,
    scalar& quadArea,
    scalar& intQuadArea
)
{
    isoDebug(Info << "Enter getQuadArea" << endl;)
    isoCutter cutter(mesh_,ap_);
    DynamicList<point> pf0(2), pf1(2);
    cutter.getFaceCutPoints(p,f,f0,pf0);
    cutter.getFaceCutPoints(p,f,f1,pf1);

    label np0(pf0.size()), np1(pf1.size());
    isoDebug(Info << "Face " << fLabel << " was cut at " << pf0 << " by f0 = " << f0 << " and at " << pf1 << " by " << " f1 = " << f1 << endl;)

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
        else if ( np0 == 1 )
        {
            B = A + 1e-4*(pf1[1]-pf1[0]);
        }
        else
        {
            Info << "Warning: Vertex face was cut at " << pf0 << " by f0 = " << f0 << endl;
        }

        if (np1 == 2 && mag(pf1[0]-pf1[1]) > SMALL)
        {
            D = pf1[1];
        }
        else if ( np1 == 1 )
        {
            D = C + 1e-4*(A-B);
        }
        else
        {
            Info << "Warning: Vertex face was cut at " << pf1 << " by f1 = " << f1 << endl;
        }

        vector zhat;
        //Defining local coordinates for area integral calculation
        if ( (mesh_.faces()[fLabel]).size() == p.size() ) //full face
        {
            zhat = mesh_.Sf()[fLabel]/mag(mesh_.Sf()[fLabel]);
        }
        else //triangle subface - make sure that triangle has same orientation as parent face
        {
            zhat = (p[1]-p[0])^(p[2]-p[0]);
            zhat /= mag(zhat);
        }

        vector xhat = B-A;
        xhat -= (xhat & zhat)*zhat;
        xhat /= mag(xhat);

        vector yhat = zhat ^ xhat;
        yhat /= mag(yhat); //Should not be necessary

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
        intQuadArea = alpha/3.0 + 0.5*beta;

////        Info << "Bx = " << Bx << ", Cx = " << Cx << ", Cy = " << Cy << ", Dx = " << Dx << ", Dy = " << Dy << ", alpha = " << alpha << ", beta = " << beta << endl;
        //area(t) = A*t^2+B*t
        //integratedArea = A/3+B/2
    }
    else
    {
        Info << "Vertex face was cut at " << pf0 << " by f0 = " << f0 << " and at " << pf1 << " by " << " f1 = " << f1 << endl;
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
    scalarField alphaNew = alpha1_ - fvc::surfaceIntegrate(dVf);
    scalar aTol = 1.0e-12;
    scalar maxAlphaMinus1 = max(alphaNew-1);
    scalar minAlpha = min(alphaNew);
    label nUndershoots = sum(neg(alphaNew+aTol));
    label nOvershoots = sum(pos(alphaNew-1-aTol));

    for ( label n = 0; n < nAlphaBounds_; n++ )
    {
        Info << "Running bounding number " << n+1 << " of time " << mesh_.time().value() << endl;

        if ( maxAlphaMinus1 > aTol )
        {
            isoDebug(Info << "Bound from above... " << endl;)
//          scalarField alpha1 = alpha1_.internalField();
//          scalarField dVfcorrected = dVf.internalField();
            surfaceScalarField dVfcorrected("dVfcorrected", dVf);
            DynamicList<label> correctedFaces(3*nOvershoots);
            boundFromAbove(alpha1_,dt,dVfcorrected,correctedFaces);
            forAll(correctedFaces,fi)
            {
                label fLabel = correctedFaces[fi];
                //Change to treat boundaries consistently
                faceValue(dVf,fLabel,faceValue(dVfcorrected,fLabel));
            }
        }

        if ( minAlpha < -aTol )
        {
            isoDebug(Info << "Bound from below... " << endl;)
            volScalarField alpha2("alpha2",1.0 - alpha1_);
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
        }

        //Check if still unbounded
        alphaNew = alpha1_ - fvc::surfaceIntegrate(dVf);
        maxAlphaMinus1 = max(alphaNew-1);
        minAlpha = min(alphaNew);
        nUndershoots = sum(neg(alphaNew+aTol));
        nOvershoots = sum(pos(alphaNew-1-aTol));
        Info << "After bounding number " << n+1 << " of time " << mesh_.time().value() << ":" << endl;
        Info << "nOvershoots = " << nOvershoots << " with max(alphaNew-1) = " << maxAlphaMinus1 << " and nUndershoots = " << nUndershoots << " with min(alphaNew) = " << minAlpha << endl;
    }
}


void Foam::isoAdvection::boundFromAbove
(
    const volScalarField& alpha1,
    const scalar dt,
    surfaceScalarField& dVf,
    DynamicList<label>& correctedFaces
)
{
    correctedFaces.clear();
    scalar aTol = 1e-12;

    forAll(alpha1,ci)
    {
        scalar Vi = mesh_.V()[ci];
        scalar alpha1New = alpha1[ci] - netFlux(dVf,ci)/Vi;
        scalar alphaOvershoot = alpha1New - 1.0;
        scalar fluidToPassOn = alphaOvershoot*Vi;
        label nFacesToPassFluidThrough = 1;

        //First try to pass surplus fluid on to neighbour cells that are not filled and to which dVf < phi*dt
        while ( alphaOvershoot > aTol && nFacesToPassFluidThrough > 0 )
        {
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
                correctedFaces.append(fLabel);
            }
            alpha1New = alpha1[ci] - netFlux(dVf,ci)/Vi;
            alphaOvershoot = alpha1New - 1.0;
            fluidToPassOn = alphaOvershoot*Vi;
            isoDebug(Info << "\nNew alpha for cell " << ci << ": " << alpha1New << endl;)
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
//  Info << "----Enter netFlux for cell " << cLabel << endl;

    //Net volume of water leaving cell cLabel. If dV is negative, cLabel, receives water from neighbours.
    scalar dV = 0.0;

    const labelList& fLabels = mesh_.cells()[cLabel];
//  const scalarField& dVfIn = dVf.internalField();

    forAll (fLabels, fi)
    {
        const label fLabel = fLabels[fi];
        const scalar dVff = faceValue(dVf,fLabel);
//      const scalar phi = faceValue(phi_,fLabel);

//      Info << "-------face " << fLabel << endl;

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
                "void isoAdvection::getDownwindFaces(...)"
            )   << "Cannot find patch for face " << fLabel
                << abort(FatalError);
        }

        const label faceI =
            mesh_.boundaryMesh()[patchI].whichFace
            (
                fLabel
            );

        ;
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
                "void isoAdvection::getDownwindFaces(...)"
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


void Foam::isoAdvection::advect
(
    const scalar dt
)
{
    isoDebug(Info << "Enter advect" << endl;)

    surfaceScalarField dVf
    (
        IOobject
        (
            "dVf",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("vol", dimVol, 0)
    );
    timeIntegratedFlux(dt, dVf);

    //Bounding
    limitFluxes(dVf,dt);
    //For each cell sum contributions from faces with pos sign for owner and neg sign for neighbour (as if it is a flux) and divide by cell volume
    alpha1_ -= fvc::surfaceIntegrate(dVf);
    alpha1_.correctBoundaryConditions();
}


void Foam::isoAdvection::getTransportedVolume
(
    const scalar dt,
    surfaceScalarField& dVf
)
{
    isoDebug(Info << "Enter advect" << endl;)

    timeIntegratedFlux(dt, dVf);

    //Bounding
    limitFluxes(dVf,dt);
}


void Foam::isoAdvection::getxSnSdotnF
(
    surfaceVectorField& xSnSdotnF
)
{
    volPointInterpolation vpi(mesh_);
    ap_ = vpi.interpolate(alpha1_);

    vectorField& xSnSdotnFi = xSnSdotnF.internalField();

    //Construction of isosurface calculator to get access to its functionality
    isoCutter cutter(mesh_,ap_);

    //Find surface cells
    DynamicList<label> surfaceCells(ceil(mesh_.nCells()/10));
    findSurfaceCells(surfaceCells);

    //Isocutting estimate of water flux from surface cellss
    interpolationCellPoint<vector> UInterp(U_);
    forAll(surfaceCells,cellI)
    {
        const label ci = surfaceCells[cellI];

        //Make list of all cell faces out of which fluid is flowing
        DynamicList<label> downwindFaces((mesh_.cells()[ci]).size());
        getDownwindFaces(ci, downwindFaces);

        //Calculate isoFace0 centre xs0, normal ns0, and velocity U0 = U(xs0)
        scalar f0(0.0), Un0(0.0);
        vector x0(vector::zero), n0(vector::zero);
        calcIsoFace(ci,x0,n0,f0,Un0,UInterp); //This one really also should give us a0 on all faces since it is calculated anyway. Do this with a cutCell structure
        Info << "calcIsoFace for cell " << ci << " gives initial surface: \nx0 = " << x0 << ", \nn0 = " << n0 << ", \nf0 = " << f0 << ", \nUn0 = " << Un0 << endl;
//        isoDebug(Info << "calcIsoFace gives initial surface: \nx0 = " << x0 << ", \nn0 = " << n0 << ", \nf0 = " << f0 << ", \nUn0 = " << Un0 << endl;)

        //Estimating time integrated water flux through each downwinding face
        forAll(downwindFaces,fi)
        {
            const label fLabel = downwindFaces[fi];
            vector nf = mesh_.Sf()[fLabel];
            nf /= mag(nf);
            //Finding average point along cutting line
            DynamicList<point> cutPoints;
            cutter.getFaceCutPoints(fLabel,f0,cutPoints);
            Info << "Cutpoints for face " << fLabel << ": " << cutPoints << endl;
            vector xs = vector::zero;
            scalar len = 0.0;

            if (cutPoints.size() > 0)
            {
                forAll(cutPoints,pi) //There will typically only be two cutting points
                {
                    xs += cutPoints[pi];
                    Info << "xs equals " << xs << endl;
                }
                xs /= cutPoints.size();
                Info << "xs equals " << xs << endl;
                for(label pi = 0; pi < cutPoints.size()-1; pi++)
                {
                    len += mag(cutPoints[pi+1]-cutPoints[pi]);
                    Info << "len equals " << len << endl;
                }
            }
            //Finding length of cutting line
            xSnSdotnFi[fLabel] = xs*(nf & n0/mag(n0))*len;
        }
    }
}


void Foam::isoAdvection::getIsoCentreAndNormal
(
    volVectorField& Ci,
    volVectorField& Si
)
{

    volPointInterpolation vpi(mesh_);
    ap_ = vpi.interpolate(alpha1_);

    //Construction of isosurface calculator to get access to its functionality
    isoCutter cutter(mesh_,ap_);

    //Find surface cells
    DynamicList<label> surfaceCells(ceil(mesh_.nCells()/10));
    findSurfaceCells(surfaceCells);

    //Isocutting estimate of water flux from surface cellss
    label maxIter(100);
    vector subCellCtr;

    forAll(surfaceCells,cellI)
    {
        const label ci = surfaceCells[cellI];
        scalar f0(0.0);
        cutter.vofCutCell(ci, alpha1_[ci], vof2IsoTol_, maxIter, f0, subCellCtr);
        vector x0(vector::zero), n0(vector::zero);
        cutter.isoFaceCentreAndArea(ci,f0,x0,n0); //Stupid to recalculate this here - should be provided by vofCutCell above

        Ci[ci] = x0;
        Si[ci] = n0;
    }
}


// ************************************************************************* //