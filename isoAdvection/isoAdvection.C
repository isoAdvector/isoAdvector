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

Foam::isoAdvection::isoAdvection
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U
)
:
    mesh_(alpha1.mesh()),
    alpha1_(alpha1),
    ap_(mesh_.nPoints(),0.0),
    phi_(phi),
    U_(U),
    isSurfaceCell_(mesh_.nCells(),false)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::isoAdvection::timeIntegratedFlux
(
    const scalar& dt,
    scalarField& dVf
)
{
    dVf = 0.0; //estimated total water volume transported across mesh faces during time interval dt

    //Interpolating VOF field to mesh points
    volPointInterpolation vpi(mesh_);
    ap_ = vpi.interpolate(alpha1_);

    //Construction of isosurface calculator to get access to its functionality
    isoCutter cutter(mesh_,ap_);

    //Find surface cells
    DynamicList<label> surfaceCells;
    findSurfaceCells(surfaceCells);

    forAll(surfaceCells,cellI)
    {
        const label ci = surfaceCells[cellI];
        Info << "\n------------ Cell " << ci << " with alpha1_ = " << alpha1_[ci] << " ------------" << endl;

        //Make list of all cell faces out of which fluid is flowing
        DynamicList<label> outFluxingFaces;
        getOutFluxFaces(ci, outFluxingFaces);
        Info << "outFluxingFaces: " << outFluxingFaces << endl;

        //Calculate isoFace0 centre xs0, normal ns0, and velocity U0 = U(xs0)
        scalar f0(0.0), Un0(0.0);
        vector x0(vector::zero), n0(vector::zero);
        calcIsoFace(ci,x0,n0,f0,Un0); //This one really also should give us a0 on all faces since it is calculated anyway. Do this with a cutCell structure

        //Estimating time integrated water flux through each outfluxing face
        forAll(outFluxingFaces,fi)
        {
            label fLabel = outFluxingFaces[fi];
            dVf[fLabel] = timeIntegratedFlux(fLabel,x0,n0,Un0,f0,dt);
//          scalar olddVf = phi_[fLabel]*0.5*(a0 + a1)*(t1-t0);
//          scalar newdVf = phi_[fLabel]/mag(mesh_.Sf()[fLabel])*(t1-t0)*(a0*mag(mesh_.Sf()[fLabel]) + sign(Un0)*integratedArea(fLabel,f0,f1));
        }
    }
    //Here go through all influx faces of surface cells and set dVf value according to alpha1 in other cell if this other cell is not a surface cell.
}

void Foam::isoAdvection::findSurfaceCells
(
    DynamicList<label>& surfaceCells
)
{
    isSurfaceCell_ = false;
    forAll(alpha1_,ci)
    {
        scalar aMin, aMax;
        subSetExtrema(ap_,mesh_.cellPoints()[ci],aMin,aMax);
        if ( (aMin < 0.5 && aMax > 0.5) && ( 1e-10 < alpha1_[ci] && alpha1_[ci] < 1-1e-10 ) )
        {
            isSurfaceCell_[ci] = true;
            surfaceCells.append(ci);
        }
    }
    surfaceCells.shrink();
}

void Foam::isoAdvection::calcIsoFace
(
    const label& ci,
    vector& x0,
    vector& n0,
    scalar& f0,
    scalar& Un0
)
{
    //Construction of isosurface calculator to get access to its functionality
    isoCutter cutter(mesh_,ap_);
    scalar tol(1e-6);
    label maxIter(100);
    vector subCellCtr;
    cutter.vofCutCell(ci, alpha1_[ci], tol, maxIter, f0, subCellCtr);

//  bool cellAlmostEmpty = mag(f0-aMax) < 1e-6*(aMax-aMin);
//  bool cellAlmostFull = mag(f0-aMin) < 1e-6*(aMax-aMin);
    bool cellAlmostEmpty = mag(alpha1_[ci]) < tol;
    bool cellAlmostFull = mag(alpha1_[ci]) > 1.0-tol;
    if ( !cellAlmostFull & !cellAlmostEmpty)
    {
        cutter.isoFaceCentreAndArea(ci,f0,x0,n0); //Stupid to recalculate this here - should be provided by vofCutCell above
        //Make n0 point out of subCell
        Info << "x0: " << x0 << ", n0: " << n0 << ", subCellCtr: " << subCellCtr << endl;
        if (((x0 - subCellCtr) & n0) < 0) //n0 should point out of water surface, i.e. in the direction from subCell centre to isoFace centre
        {
            n0 *= (-1.0);
            Info << "Changing direction of n0 " << n0 << endl;
        }
    }
    else if (cellAlmostEmpty) //We need to get the isoFace normal from an isoFace which is just inside the cell
    {
        scalar f0Inside =  f0 - tol;
        cutter.isoFaceCentreAndArea(ci,f0Inside,x0,n0);
        if (((mesh_.C()[ci] - x0) & n0) < 0.0) //n0 should point from isoFace centre towards cell centre for an almost empty cell
        {
            n0 *= (-1.0);
            Info << "Changing direction of n0 " << n0 << endl;
        }
    }
    else //cellAlmostFull
    {
        scalar f0Inside = f0 + tol;
        cutter.isoFaceCentreAndArea(ci,f0Inside,x0,n0);
        if (((x0 - mesh_.C()[ci]) & n0) < 0.0) //n0 should point from cell centre towards isoFace centre for an almost full cell
        {
            n0 *= (-1.0);
            Info << "Changing direction of n0 " << n0 << endl;
        }
    }
    if (mag(n0) > 1e-10)
    {
        n0 /= mag(n0);
    }
    else
    {
        Info << "Warning: mag(n0) = " << mag(n0) << ". Cannot normalise it." << endl;
    }

    //Interpolate velocity to isoFace centre
    interpolationCellPoint<vector> UInterp(U_);
    vector U0 = UInterp.interpolate(x0,ci);
    Info << "U0 = " << U0 << endl;
    Un0 = (U0 & n0);
    Info << "Un0 = " << Un0 << endl;
    Info << "Initial values for time step:" << endl;
    Info << "f0 = " << f0 << ", subCellCentre = " << subCellCtr << ", isoFaceCentre = " << x0 << ", isoFaceNormal = " << n0 << ", isoFaceVelocity = " << U0 << ", isoFaceNormalVelocity = " << Un0 << endl;

}

Foam::scalar Foam::isoAdvection::timeIntegratedFlux
(
    const label& fLabel,
    const vector& x0,
    const vector& n0,
    const scalar& Un0,
    const scalar& f0,
    const scalar& dt
)
{
    //Find sorted list of times where the isoFace will arrive at face points given initial position x0 and velocity Un0*n0
    labelList pLabels = mesh_.faces()[fLabel];
    label nPoints = pLabels.size();
    pointField fPts(nPoints);
    forAll(fPts,pi)
    {
        fPts[pi] = mesh_.points()[pLabels[pi]];
    }
    scalarField pTimes = ((fPts - x0) & n0)/Un0; //Here we estimate time of arrival to the face points from their normal distance to the initial surface and the surface normal velocity
    scalarField sortedTimes(pTimes);
    sort(sortedTimes);

    scalar dVf(0.0); //Volume flowing through face in time interval [0,dt] to be calculated below

    //Dealing with case where face is not cut by surface during time interval [0,dt] because face was already passed by surface
    if ( sortedTimes[nPoints-1] <= 0.0 ) //All cuttings in the past
    {
        dVf = phi_[fLabel]*dt*sign(Un0); //If all face cuttings were in the past and cell is filling up (Un0>0) then face must be full during whole time interval
        return dVf;
    }

    //Dealing with case where face is not cut by surface during time interval [0,dt] because dt is too small for surface to reach closest face point
    if ( sortedTimes[0] >= dt ) //All cuttings in the future
    {
        dVf = phi_[fLabel]*dt*(1-sign(Un0)); //If all cuttings are in the future but non of them within [0,dt] then if cell is filling up (Un0 > 0) face must be empty during whole time interval
        return dVf;
    }

    //Cutting sortedTimes at 0 and dt
    scalarField t;
    t.append(0.0);
    forAll(sortedTimes,ti)
    {
        if ( 0.0 < sortedTimes[ti] && sortedTimes[ti] < dt )
        {
            t.append(sortedTimes[ti]);
        }
    }
    t.append(dt);

    bool faceUncutInFirstInterval(sortedTimes[0] > 0.0);
    bool faceUncutInLastInterval(sortedTimes[nPoints-1] < dt);

    if (t.size() < 3)
    {
        Info << "Warning for face " << fLabel << ": tIntevals contains less that 3 points!" << endl;
    }

    //Dealing with cases where face is cut at least once during time interval [0,dt]
    label nt = 0; //Time interval counter
    scalar initialArea(0.0); //Submerged area at time
    if ( faceUncutInFirstInterval ) //Special treatment for first time interval if face is uncut during this
    {
        dVf = phi_[fLabel]*(t[nt+1]-t[nt])*(1.0-sign(Un0)); //If Un0 > 0 cell is filling up - hence if face is cut at a later time but not initially it must be initially empty
        initialArea = mag(mesh_.Sf()[fLabel])*(1.0-sign(Un0));
        nt++;
    }

    while ( nt < t.size()-(1+faceUncutInLastInterval) ) //
    {
        scalar A(0.0), B(0.0);
        quadAreaCoeffs(fLabel,pTimes,t[nt],t[nt+1],A,B);
        scalar integratedQuadArea = sign(Un0)*(A/3.0 + 0.5*B); //Integration of area(t) = A*t^2+B*t from t = 0 to 1.
        scalar Unf = phi_[fLabel]/mag(mesh_.Sf()[fLabel]); //normal velocity on face
        dVf += Unf*(t[nt+1]-t[nt])*(initialArea + integratedQuadArea);
        initialArea += sign(Un0)*(A + B); //Adding quad area to
//      scalar newdVf = phi_[fLabel]/mag(mesh_.Sf()[fLabel])*(t1-t0)*(a0*mag(mesh_.Sf()[fLabel]) + sign(Un0)*integratedArea(fLabel,f0,f1));
        nt++;
    }

    if ( faceUncutInLastInterval ) //Special treatment for last time interval if face is uncut during this
    {
        dVf += phi_[fLabel]*(t[nt+1]-t[nt])*sign(Un0); //If face is cut at some intermediate time but not at last time, then if Un0 > 0 (cell filling up) face must be filled at last time interval.
    }

    return dVf;
}


void Foam::isoAdvection::getOutFluxFaces
(
    const label& ci,
    DynamicList<label>& outFluxingFaces
)
{
    const cellList& cells = mesh_.cells();
    forAll(cells[ci],fi)
    {
        label fLabel = cells[ci][fi];
        if (mesh_.owner()[fLabel] == ci)
        {
            if (phi_[fLabel] > 1e-10) //HARDCODED LIMIT
            {
                outFluxingFaces.append(fLabel);
            }
        }
        else if ( phi_[fLabel] < -1e-10 ) //ci must be neighbour of fLabel
        {
            outFluxingFaces.append(fLabel);
        }
    }
    outFluxingFaces.shrink();
}


Foam::scalar Foam::isoAdvection::integratedArea
(
    const label& fLabel,
    const scalar& f0,
    const scalar& f1
)
{
    isoCutter cutter(mesh_,ap_);
    pointField pf0, pf1;
    cutter.getFaceCutPoints(fLabel,f0,pf0);
    cutter.getFaceCutPoints(fLabel,f1,pf1);
    label np0(pf0.size()), np1(pf1.size());
    Info << "Face " << fLabel << " was cut " << np0 << " times by f0 = " << f0 << " and " << np1 << " times by " << " f1 = " << f1 << endl;

    scalar area(0.0);

    if ( np0 > 0 && np1 > 0)
    {
        //Defining local coordinates for area integral calculation
        vector xhat(vector::zero), yhat, zhat;
        zhat = mesh_.Sf()[fLabel]/mag(mesh_.Sf()[fLabel]);
        if (np0 == 2)
        {
            xhat = pf0[1]-pf0[0];
            xhat -= (xhat & zhat)*zhat;
            xhat /= mag(xhat);
        }
        else if (np1 == 2)
        {
            xhat = pf1[1]-pf1[0];
            xhat -= (xhat & zhat)*zhat;
            xhat /= mag(xhat);
        }
        yhat = zhat ^ xhat;
        yhat /= mag(yhat); //Should not be necessary

        Info << "xhat = " << xhat << ", yhat = " << yhat << ", zhat = " << zhat << endl;//". x.x = " << xhat & xhat << ", y.y = " << yhat & yhat <<", z.z = " << zhat & zhat << ", x.y = " << xhat & yhat << ", x.z = " << xhat & zhat << ", y.z = " << yhat & zhat << endl;

        //Triangle cases
        if (np0 == 1)
        {
            pf0.append(pf0[0]);
        }
        else if (np1 == 1)
        {
            pf1.append(pf1[0]);
        }

        //Swapping pf1 points if pf0 and pf1 point in same general direction (because we want a quadrilateral ABCD where pf0 = AB and pf1 = CD)
        if ( ((pf0[1]-pf0[0]) & (pf1[1]-pf1[0])) > 0 )
        {
            vector tmp = pf1[0];
            pf1[0] = pf1[1];
            pf1[1] = tmp;
        }

        scalar Bx = mag(pf0[1]-pf0[0]);
        scalar Cx = (pf1[0]-pf0[0]) & xhat;
        scalar Cy = mag((pf1[0]-pf0[0]) & yhat);
        scalar Dx = (pf1[1]-pf0[0]) & xhat;
        scalar Dy = mag((pf1[1]-pf0[0]) & yhat);

//      Info << "Bx = " << Bx << ", Cx = " << Cx << ", Cy = " << Cy << ", Dx = " << Dx << ", Dy = " << Dy << endl;

        area = ((Cx-Bx)*Dy-Dx*Cy)/6.0 + 0.25*Bx*(Dy+Cy);
    }

    return area;
}


void Foam::isoAdvection::quadAreaCoeffs
(
    const label& fLabel,
    const scalarField& f,
    const scalar& f0,
    const scalar& f1,
    scalar& A,
    scalar& B
)
{
    isoCutter cutter(mesh_,ap_);
    pointField pf0, pf1;
    cutter.getFaceCutPoints(fLabel,f,f0,pf0);
    cutter.getFaceCutPoints(fLabel,f,f1,pf1);
//  cutter.getFaceCutPoints(fLabel,f0,pf0);
//  cutter.getFaceCutPoints(fLabel,f1,pf1);

    label np0(pf0.size()), np1(pf1.size());
    Info << "Face " << fLabel << " was cut " << np0 << " times by f0 = " << f0 << " and " << np1 << " times by " << " f1 = " << f1 << endl;

//  scalar area(0.0);
    A = 0.0;
    B = 0.0;

    if ( np0 > 0 && np1 > 0)
    {
        //Defining local coordinates for area integral calculation
        vector xhat(vector::zero), yhat, zhat;
        zhat = mesh_.Sf()[fLabel]/mag(mesh_.Sf()[fLabel]);
        if (np0 == 2)
        {
            xhat = pf0[1]-pf0[0];
            xhat -= (xhat & zhat)*zhat;
            xhat /= mag(xhat);
        }
        else if (np1 == 2)
        {
            xhat = pf1[1]-pf1[0];
            xhat -= (xhat & zhat)*zhat;
            xhat /= mag(xhat);
        }
        yhat = zhat ^ xhat;
        yhat /= mag(yhat); //Should not be necessary

        Info << "xhat = " << xhat << ", yhat = " << yhat << ", zhat = " << zhat << endl;//". x.x = " << xhat & xhat << ", y.y = " << yhat & yhat <<", z.z = " << zhat & zhat << ", x.y = " << xhat & yhat << ", x.z = " << xhat & zhat << ", y.z = " << yhat & zhat << endl;

        //Triangle cases
        if (np0 == 1)
        {
            pf0.append(pf0[0]);
        }
        else if (np1 == 1)
        {
            pf1.append(pf1[0]);
        }

        //Swapping pf1 points if pf0 and pf1 point in same general direction (because we want a quadrilateral ABCD where pf0 = AB and pf1 = CD)
        if ( ((pf0[1]-pf0[0]) & (pf1[1]-pf1[0])) > 0 )
        {
            vector tmp = pf1[0];
            pf1[0] = pf1[1];
            pf1[1] = tmp;
        }

        scalar Bx = mag(pf0[1]-pf0[0]);
        scalar Cx = (pf1[0]-pf0[0]) & xhat;
        scalar Cy = mag((pf1[0]-pf0[0]) & yhat);
        scalar Dx = (pf1[1]-pf0[0]) & xhat;
        scalar Dy = mag((pf1[1]-pf0[0]) & yhat);

//      Info << "Bx = " << Bx << ", Cx = " << Cx << ", Cy = " << Cy << ", Dx = " << Dx << ", Dy = " << Dy << endl;

//      area = ((Cx-Bx)*Dy-Dx*Cy)/6.0 + 0.25*Bx*(Dy+Cy);
        A = 0.5*((Cx-Bx)*Dy-Dx*Cy);
        B = 0.5*Bx*(Dy+Cy);
        //area(t) = A*t^2+B*t
        //integratedArea = A/3+B/2
    }
//  return area;
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

void Foam::isoAdvection::boundAlpha
(
    scalarField& dVf
)
{
    const scalarField& V = mesh_.V();
    scalarField dV(alpha1_.size(),0.0);

    const scalar tol = 1e-8;//SMALL;
    Info << "Bounding transport (tolerance = " << tol << ")" << endl;

    DynamicList<label> overFilledCells, overEmptiedCells;
    forAll(alpha1_,ci)
    {
        dV[ci] = netFlux(dVf,ci);
        bool cellIsOverFilled = dV[ci] < -tol && mag(dV[ci]) > V[ci]*(1.0-alpha1_[ci]) + tol;
        if ( cellIsOverFilled )
        {
//          Info << "Cell " << ci << " is overfilled: -dV = " << -dV[ci] << " < V*(1-alpha) = " << V[ci]*(1.0-alpha1_[ci]) << endl;
            overFilledCells.append(ci);
        }
        else
        {
            bool cellIsOverEmptied = dV[ci] > V[ci]*alpha1_[ci] + tol;
            if ( cellIsOverEmptied )
            {
//              Info << "Cell " << ci << " is overemptied: dV = " << dV[ci] << " > V*alpha = " << V[ci]*alpha1_[ci] << endl;
                overEmptiedCells.append(ci);
            }
        }
    }
    overFilledCells.shrink();
    overEmptiedCells.shrink();

    //Correcting overfilling fluxes
    {
        label nOverFilledCells = overFilledCells.size();
        label ci = 0;
        while(ci < nOverFilledCells)
        {
            const label cLabel = overFilledCells[ci];
            DynamicList<label> outFluxingFaces;
            getOutFluxFaces(cLabel, outFluxingFaces);
            scalar sumdVOut = 0.0;
            forAll(outFluxingFaces,fi) //If some outfluxing faces e.g. those fluxing pure air or water are to be left out, this is the place to do it.
            {
                const label fLabel = outFluxingFaces[fi];
                sumdVOut += mag(dVf[fLabel]);
            }
    //      scalar surplusWater = V[cLabel]*(alpha1_[cLabel]-1.0) - dV[cLabel];
            scalar surplusWater = mag(dV[cLabel])-V[cLabel]*(1.0-alpha1_[cLabel]);
            const scalar alphaNew = alpha1_[cLabel]-dV[cLabel]/V[cLabel];
            Info << "\nCell " << cLabel << " has alpha = " << alpha1_[cLabel] << ", dV = " << dV[cLabel] << " m3, surplusWater = " << surplusWater << " m3 which would give alphaNew = " << alphaNew << endl;
            if (surplusWater>tol)
            {
                Info << "sumdVOut = " << sumdVOut << endl;
                //          Info << "\nCell " << cLabel << " has surplusWater = " << surplusWater << " m3." << endl;
                forAll(outFluxingFaces,fi)
                {
                    const label fLabel = outFluxingFaces[fi];
                    scalar waterAdded = surplusWater*dVf[fLabel]/sumdVOut;
                    dVf[fLabel] += waterAdded; //Whatever sign dVf has, it corresponds to out of cLabel. We must make it smaller, thus adding a little if it is negative and subtracting a little if it is positive
                    //Checking if receiving cell becomes overfilled by the additional water it received
                    label otherCellLabel;
                    if (cLabel == mesh_.owner()[fLabel])
                    {
                        otherCellLabel = mesh_.neighbour()[fLabel];
                    }
                    else
                    {
                        otherCellLabel = mesh_.owner()[fLabel];
                    }
                    scalar newdVOther = dV[otherCellLabel] - mag(waterAdded);
                    bool otherCellIsOverfilled = newdVOther < -tol && mag(newdVOther) > V[otherCellLabel]*(1.0-alpha1_[otherCellLabel]) + tol;

                    if ( otherCellIsOverfilled )
                    {
                        bool otherCellWasAlreadyOverfilled = dV[otherCellLabel] < -tol && (mag(dV[otherCellLabel]) > V[otherCellLabel]*(1.0-alpha1_[otherCellLabel]) + tol);
                        if ( !otherCellWasAlreadyOverfilled )
                        {
                            Info << "Cell " << otherCellLabel << " was overfilled. Appending to overFilledCells list." << endl;
                            overFilledCells.append(otherCellLabel);
                            nOverFilledCells++;
                        }
                    }
                    dV[otherCellLabel] = newdVOther; //Updating dV for neighbour cell - Do not move this in front of the above if loop
    //              dV[otherCellLabel] = netFlux(dVf,otherCellLabel);
                    Info << "Moving " << waterAdded << "  m3 to its neighbour cell " << otherCellLabel << "." << endl;
                }
    //          dV[cLabel] = netFlux(dVf,cLabel);
                dV[cLabel] += surplusWater; //Updating dV for overfilled cell
            }
            ci++;
        }
    }
////////////////////
    //Correcting overemptying fluxes
    {
        label nOverEmptiedCells = overEmptiedCells.size();
        label ci = 0;
        while(ci < nOverEmptiedCells)
        {
            const label cLabel = overEmptiedCells[ci];
            DynamicList<label> outFluxingFaces;
            getOutFluxFaces(cLabel, outFluxingFaces);
            scalar sumdVOut = 0.0;
            forAll(outFluxingFaces,fi) //If some outfluxing faces e.g. those fluxing pure air or water are to be left out, this is the place to do it.
            {
                const label fLabel = outFluxingFaces[fi];
                sumdVOut += mag(dVf[fLabel]);
            }
            scalar missingWater = dV[cLabel] - V[cLabel]*alpha1_[cLabel];
            const scalar alphaNew = alpha1_[cLabel]-dV[cLabel]/V[cLabel];
            Info << "\nCell " << cLabel << " has alpha = " << alpha1_[cLabel] << ", dV = " << dV[cLabel] << " m3, missingWater = " << missingWater << " m3 which would give alphaNew = " << alphaNew << endl;
            if (missingWater>tol)
            {
                Info << "sumdVOut = " << sumdVOut << endl;
                forAll(outFluxingFaces,fi)
                {
                    const label fLabel = outFluxingFaces[fi];
                    scalar waterSubtracted = missingWater*dVf[fLabel]/sumdVOut;
                    dVf[fLabel] -= waterSubtracted; //Whatever sign dVf has, it corresponds to out of cLabel. We must make it smaller, thus adding a little if it is negative and subtracting a little if it is positive
                    //Checking if receiving cell becomes overemptied by the decrease in water it receives
                    label otherCellLabel;
                    if (cLabel == mesh_.owner()[fLabel])
                    {
                        otherCellLabel = mesh_.neighbour()[fLabel];
                    }
                    else
                    {
                        otherCellLabel = mesh_.owner()[fLabel];
                    }
                    scalar newdVOther = dV[otherCellLabel] + mag(waterSubtracted);
                    bool otherCellIsOverEmptied = newdVOther > V[otherCellLabel]*alpha1_[otherCellLabel] + tol;
                    if ( otherCellIsOverEmptied )
                    {
                        bool otherCellWasAlreadyOverEmptied = dV[otherCellLabel] > V[otherCellLabel]*alpha1_[otherCellLabel] + tol;
                        if ( !otherCellWasAlreadyOverEmptied )
                        {
                            Info << "Cell " << otherCellLabel << " was overemptied. Appending to overEmptiedCells list." << endl;
                            overEmptiedCells.append(otherCellLabel);
                            nOverEmptiedCells++;
                        }
                    }
                    dV[otherCellLabel] = newdVOther; //Updating dV for neighbour cell - Do not move this in front of the above if loop
    //              dV[otherCellLabel] = netFlux(dVf,otherCellLabel);
                    Info << "Moving " << waterSubtracted << "  m3 from its neighbour cell " << otherCellLabel << "." << endl;
                }
                dV[cLabel] -= missingWater; //Updating dV for overemptied cell
    //          dV[cLabel] = netFlux(dVf,cLabel);
            }
            ci++;
        }
    }
}

Foam::scalar Foam::isoAdvection::netFlux
(
    const scalarField& dVf,
    const label& cLabel
)
{
    scalar dV = 0.0;
    const labelList& fLabels = mesh_.cells()[cLabel];
    forAll(fLabels,fi)
    {
        if (fLabels[fi] < mesh_.nInternalFaces())  //Error when this is removed - should be removed when boundary face treatment is implemented
        {
            if (mesh_.owner()[fLabels[fi]] == cLabel)
            {
                dV += dVf[fLabels[fi]];
            }
            else
            {
                dV -= dVf[fLabels[fi]];
            }
        }
    }
    return dV;
}

void Foam::isoAdvection::advect
(
    const scalar& dt
)
{
    surfaceScalarField dVf(0*phi_); //Construct as copy - has same dimensions so will probably give problems. How to construct as zero's with phi's mesh etc e.g. surfaceScalarField dVf(phi.size(),0.0)?
    dVf.dimensions().reset(phi_.mesh().V().dimensions());
    scalarField& dVfi = dVf;
    timeIntegratedFlux(dt, dVfi);
//  boundAlpha(dVfi);
    alpha1_ -= fvc::surfaceIntegrate(dVf); //For each cell sum contributions from faces with pos sign for owner and neg sign for neighbour (as if it is a flux) and divide by cell volume
    alpha1_.correctBoundaryConditions();
//  alpha1_ = min(1.0,max(0.0,alpha1_));
}


// ************************************************************************* //
