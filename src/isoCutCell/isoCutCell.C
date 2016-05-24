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

#include "isoCutCell.H"

#ifdef DETAILS2LOG
#define isoDebug(x) x
#else
#define isoDebug(x)
#endif

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isoCutCell::isoCutCell
(
    const fvMesh& mesh,
    scalarField& f
)
:
    mesh_(mesh),
    cellI_(-1),
    f_(f),
    isoValue_(0),
    isoCutFace_(isoCutFace(mesh_,f_)),
    isoCutFaces_(10),
    isoCutFacePoints_(10),
    isoCutFaceCentres_(10),
    isoCutFaceAreas_(10),
    isoFaceEdges_(10),
    isoFacePoints_(10),
    isoFaceCentre_(vector::zero),
    isoFaceArea_(vector::zero),
    subCellCentre_(vector::zero),
    subCellVolume_(-10),
    VOF_(-10),
    fullySubFaces_(10),
    cellStatus_(-1),
    subCellCentreAndVolumeCalculated_(false),
    isoFaceCentreAndAreaCalculated_(false)
{
    clearStorage();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::isoCutCell::calcSubCellCentreAndVolume()
{
    if (cellStatus_ == 0)
    {
        // Clear the fields for accumulation
        subCellCentre_ = vector::zero;
        subCellVolume_ = 0.0;

        // first estimate the approximate cell centre as the average of
        // face centres

        vector cEst = isoFaceCentre_;

        forAll(isoCutFaceCentres_, facei)
        {
            cEst += isoCutFaceCentres_[facei];
        }

        forAll(fullySubFaces_, facei)
        {
            cEst += mesh_.Cf()[fullySubFaces_[facei]];
        }

        label nCellFaces(1 + isoCutFaceCentres_.size() + fullySubFaces_.size());
        cEst /= nCellFaces;

        //Contribution to subcell centre and volume from isoface
        scalar pyr3Vol =
            max(mag(isoFaceArea_ & (isoFaceCentre_ - cEst)), VSMALL);

        // Calculate face-pyramid centre
        vector pc = (3.0/4.0)*isoFaceCentre_ + (1.0/4.0)*cEst;

        // Accumulate volume-weighted face-pyramid centre
        subCellCentre_ += pyr3Vol*pc;

        // Accumulate face-pyramid volume
        subCellVolume_ += pyr3Vol;
        
        //Contribution to subcell centre and volume from cut faces
        forAll(isoCutFaceCentres_, facei)
        {
            // Calculate 3*face-pyramid volume
            scalar pyr3Vol =
                max(mag(isoCutFaceAreas_[facei] & (isoCutFaceCentres_[facei] - cEst)), VSMALL);

            // Calculate face-pyramid centre
            vector pc = (3.0/4.0)*isoCutFaceCentres_[facei] + (1.0/4.0)*cEst;

            // Accumulate volume-weighted face-pyramid centre
            subCellCentre_ += pyr3Vol*pc;

            // Accumulate face-pyramid volume
            subCellVolume_ += pyr3Vol;
        }

        //Contribution to subcell centre and volume from fully submerged faces
        forAll(fullySubFaces_, facei)
        {
            const point& fCentre = mesh_.Cf()[fullySubFaces_[facei]]; 
            const vector& fArea = mesh_.Sf()[fullySubFaces_[facei]]; 

            // Calculate 3*face-pyramid volume
            scalar pyr3Vol =
                max(mag(fArea & (fCentre - cEst)), VSMALL);

            // Calculate face-pyramid centre
            vector pc = (3.0/4.0)*fCentre + (1.0/4.0)*cEst;

            // Accumulate volume-weighted face-pyramid centre
            subCellCentre_ += pyr3Vol*pc;

            // Accumulate face-pyramid volume
            subCellVolume_ += pyr3Vol;
        }

        subCellCentre_ /= subCellVolume_;
        subCellVolume_ *= (1.0/3.0);
        VOF_ = subCellVolume_/mesh_.V()[cellI_];
        
        subCellCentreAndVolumeCalculated_ = true;
/*        
        //Testing
        vector sumSf = isoFaceArea_;
        scalar sumMagSf = mag(isoFaceArea_);
        forAll(isoCutFaceCentres_, facei)
        {
            sumSf += isoCutFaceAreas_[facei];
            sumMagSf += mag(isoCutFaceAreas_[facei]);
        }
        forAll(fullySubFaces_, facei)
        {
            sumSf += mesh_.Sf()[fullySubFaces_[facei]]; 
            sumMagSf += mag(isoCutFaceAreas_[facei]);
        }
        if (mag(sumSf) > 1e-10)
        {
            Info << "Warninig: mag(sumSf)/magSumSf = " << mag(sumSf)/sumMagSf << " for surface cell" 
                << cellI_ << endl;
        }
*/
    }
    else if (cellStatus_ == 1) //cell fully above isosurface
    {
        subCellCentre_ = vector::zero;
        subCellVolume_ = 0;
        VOF_ = 0;       
    }
    else if (cellStatus_ == -1) //cell fully below isosurface
    {
        subCellCentre_ = mesh_.C()[cellI_];
        subCellVolume_ = mesh_.V()[cellI_];
        VOF_ = 1;         
    }
}


void Foam::isoCutCell::calcIsoFaceCentreAndArea()
{
    //Initial guess of face centre from edge points
    point fCentre = vector::zero;
    label nEdgePoints = 0;
    forAll(isoFaceEdges_, ei)
    {
        DynamicList<point>& edgePoints = isoFaceEdges_[ei];
        forAll(edgePoints, pi)
        {
            fCentre += edgePoints[pi];
            nEdgePoints++;
        }        
    }
    if (nEdgePoints > 0)
    {
        fCentre /= nEdgePoints;
    }
    else
    {
        Info << "Warning: nEdgePoints = 0 for cell " << cellI_ << endl;
    }
    
    vector sumN = vector::zero;
    scalar sumA = 0.0;
    vector sumAc = vector::zero;

    forAll(isoFaceEdges_, ei)
    {
        const DynamicList<point>& edgePoints = isoFaceEdges_[ei];
        const label nPoints = edgePoints.size();
        for (label pi = 0; pi < nPoints-1; pi++)
        {
            const point& nextPoint = edgePoints[pi + 1];

            vector c = edgePoints[pi] + nextPoint + fCentre;
            vector n = (nextPoint - edgePoints[pi])^(fCentre - edgePoints[pi]);
            scalar a = mag(n);

            sumN += Foam::sign(n & sumN)*n; //Edges may have different orientation
            sumA += a;
            sumAc += a*c;
        }        
    }
    
    // This is to deal with zero-area faces. Mark very small faces
    // to be detected in e.g., processorPolyPatch.
    if (sumA < ROOTVSMALL)
    {
        isoFaceCentre_ = fCentre;
        isoFaceArea_ = vector::zero;
    }
    else
    {
        isoFaceCentre_ = (1.0/3.0)*sumAc/sumA;
        isoFaceArea_ = 0.5*sumN;
    }

    
    //Check isoFaceArea_ direction and change if not pointing out of subcell
    if ( (isoFaceArea_ & (isoFaceCentre_ - subCellCentre())) < 0 )
    {
        isoFaceArea_ *= (-1);
    }
    
    isoFaceCentreAndAreaCalculated_ = true;
}


void Foam::isoCutCell::calcIsoFacePointsFromEdges()
{
//    Info << "Enter calcIsoFacePointsFromEdges() with isoFaceArea_ = " 
//        << isoFaceArea_ << " and isoFaceCentre_ = " << isoFaceCentre_ 
//        << " and isoFaceEdges_ = " << isoFaceEdges_ << endl;
    //Defining local coordinates with zhat along isoface normal and xhat from 
    //isoface centre to first point in isoFaceEdges_
    const vector zhat = isoFaceArea_/mag(isoFaceArea_);
    vector xhat = isoFaceEdges_[0][0]-isoFaceCentre_;
    xhat = (xhat - (xhat & zhat)*zhat);
    xhat /= mag(xhat);
    vector yhat = zhat^xhat;
    yhat /= mag(yhat);

//    Info << "Calculated local coordinates" << endl;
    
    //Calculating isoface point angles in local coordinates
    DynamicList<point> unsortedIsoFacePoints(3*isoFaceEdges_.size());
    DynamicList<scalar> unsortedIsoFacePointAngles(3*isoFaceEdges_.size());
    forAll(isoFaceEdges_, ei)
    {
        const DynamicList<point>& edgePoints = isoFaceEdges_[ei];
        forAll(edgePoints,pi)
        {
            const point p = edgePoints[pi];
            unsortedIsoFacePoints.append(p);
            unsortedIsoFacePointAngles.append
            (
                Foam::atan2
                (
                    ((p-isoFaceCentre_) & yhat), 
                    ((p-isoFaceCentre_) & xhat)
                )
            );
        }
    }
    
//    Info << "Calculated isoFace point angles" << endl;

    //Sorting isoface points by angle and inserting into isoFacePoints_
    labelList order(unsortedIsoFacePointAngles.size());
    Foam::sortedOrder(unsortedIsoFacePointAngles, order);
    isoFacePoints_.append(unsortedIsoFacePoints[order[0]]);
    for(label pi = 1; pi < order.size(); pi++)
    {
        if 
        (
            mag(unsortedIsoFacePointAngles[order[pi]]
                -unsortedIsoFacePointAngles[order[pi-1]]) > 1e-8
        )
        {
            isoFacePoints_.append(unsortedIsoFacePoints[order[pi]]);
        }
    }
//    Info << "Sorted isoface points by angle" << endl;

}


Foam::label Foam::isoCutCell::calcSubCell
(
    const label cellI,
    const scalar isoValue
)
{        
    //Populate isoCutFaces_, isoCutFacePoints_, fullySubFaces_, isoFaceCentres_
    //and isoFaceArea_.

    clearStorage();
    cellI_ = cellI;
    isoValue_ = isoValue;
    const labelList& faces = mesh_.cells()[cellI];
    
    forAll(faces,fi)
    {
        const label faceStatus = isoCutFace_.calcSubFace(faces[fi],isoValue_);
        if (faceStatus == 0) //face is cut
        {
            isoCutFacePoints_.append(isoCutFace_.subFacePoints());
            isoCutFaceCentres_.append(isoCutFace_.subFaceCentre());
            isoCutFaceAreas_.append(isoCutFace_.subFaceArea());
            isoFaceEdges_.append(isoCutFace_.surfacePoints());
        }
        else if (faceStatus == -1) //face fully below
        {
            fullySubFaces_.append(faces[fi]);
        }            
    }
    
    if (isoCutFacePoints_.size() > 0) //cell cut at least at one face
    {
        cellStatus_ = 0;
        calcIsoFaceCentreAndArea();
    }
    else if (fullySubFaces_.size() == 0) //cell fully above isosurface
    {
        cellStatus_ = 1;
    }
    else //cell fully below isosurface
    {
        cellStatus_ = -1;
    }
    //else all faces are below isosurface and cellStatus_ = -1,
    //which is its default value, so no action required here

    return cellStatus_;
}


Foam::point Foam::isoCutCell::subCellCentre()
{
    if ( !subCellCentreAndVolumeCalculated_ )
    {
        calcSubCellCentreAndVolume();
    }

    return subCellCentre_;
}


Foam::scalar Foam::isoCutCell::subCellVolume()
{
    if ( !subCellCentreAndVolumeCalculated_ )
    {
        calcSubCellCentreAndVolume();
    }

    return subCellVolume_;
}


Foam::DynamicList<Foam::point> Foam::isoCutCell::isoFacePoints()
{
    if (cellStatus_ == 0 && isoFacePoints_.size() == 0)
    {
        calcIsoFacePointsFromEdges();
    }
    return isoFacePoints_;
}


Foam::point Foam::isoCutCell::isoFaceCentre()
{
    if (!isoFaceCentreAndAreaCalculated_)
    {
        calcIsoFaceCentreAndArea();
    }
    return isoFaceCentre_;
}

Foam::vector Foam::isoCutCell::isoFaceArea()
{   
    if (!isoFaceCentreAndAreaCalculated_)
    {
        calcIsoFaceCentreAndArea();
    }
    return isoFaceArea_;
}
        

Foam::scalar Foam::isoCutCell::VolumeOfFluid()
{
    if ( !subCellCentreAndVolumeCalculated_ )
    {
        calcSubCellCentreAndVolume();
    }
    
    return VOF_;
}

Foam::scalar Foam::isoCutCell::isoValue()
{    
    return isoValue_;
}


void Foam::isoCutCell::clearStorage()
{
    cellI_ = -1;
    isoValue_ = 0;
    isoCutFace_.clearStorage();
    isoCutFaces_.clear();
    isoCutFacePoints_.clear();
    isoCutFaceCentres_.clear();
    isoCutFaceAreas_.clear();
    isoFaceEdges_.clear();
    isoFacePoints_.clear();
    isoFaceCentre_ = vector::zero;
    isoFaceArea_ = vector::zero;
    subCellCentre_ = vector::zero;
    subCellVolume_ = -10;
    VOF_ = -10;
    fullySubFaces_.clear();
    cellStatus_ = -1;
    subCellCentreAndVolumeCalculated_ = false;
    isoFaceCentreAndAreaCalculated_ = false;
}


Foam::scalar Foam::isoCutCell::vofCutCell
(
    const label cellI,
    const scalar alpha1,
    const scalar tol,
    const label maxIter
)
{
    //Finding cell vertex extremum values
    scalar fMin(GREAT), fMax(-GREAT);
    const labelList& pLabels = mesh_.cellPoints(cellI);
    forAll(pLabels,pi)
    {
        scalar fp = f_[pLabels[pi]];
        if (fp > fMax)
        {
            fMax = fp;
        }
        if (fp < fMin)
        {
            fMin = fp;
        }
    }

    //Handling special case where method is handed an almost full or empty cell
    if ( alpha1 < tol )
    {
        return fMax;
    }
    else if (1-alpha1 < tol)
    {
        return fMin;
    }

    //Initial guess of isovalue
    scalar aMin(0), aMax(1);
    scalar f0 = (alpha1 - aMin)/(aMax-aMin)*(fMin-fMax) + fMax;
//    scalar f0 = 0.5*(fMin + fMax);

    calcSubCell(cellI,f0);
    scalar alpha0 = VolumeOfFluid();

    //Bisection method to find root.
    //Since function is monotonically increasing and only piecewise smooth 
    //derivative based root finding algorithms should not be used here.
    label nIter = 0;
//    Info << "Cell " << cellI << ", target: alpha1 =  " << alpha1 << ":" << endl;
//    Info << "(f0 alpha) = " << f0 << " " << alpha0 << endl;
    while ( mag(alpha1-alpha0) > tol && nIter < maxIter )
    {
        if ( alpha0 > alpha1 )
        {
            aMax = alpha0;
            fMin = f0;
        }
        else
        {
            aMin = alpha0;
            fMax = f0;
        }
        f0 = 0.5*(fMin + fMax);
//        f0 = (alpha1 - aMin)/(aMax-aMin)*(fMin-fMax) + fMax;
        calcSubCell(cellI,f0);
        alpha0 = VolumeOfFluid();
//      Info << "Cell " << cellI << ", iter " << nIter << ": f0 = " << f0 << " (" << 1-f0 << ") gives alpha = " << alpha0 << endl;
//        Info << f0 << " " << alpha0 << endl;
        nIter++;
    }
//    Info << nIter-1 << ": f0 = " << f0 << " gives alpha = " << alpha0 << endl;
        
    return f0;
}


Foam::label Foam::isoCutCell::vofCutCell2
(
    const label cellI,
    const scalar alpha1,
    const scalar tol,
    const label maxIter
)
{
////    Info << "\n------ vofCutCell2 for cell " << cellI << " with alpha1 = " 
////        << alpha1 << " ------" << endl;
    //Finding cell vertex extremum values
    const labelList& pLabels = mesh_.cellPoints(cellI);
    scalarField fvert(pLabels.size());
    forAll(pLabels,pi)
    {
        fvert[pi] = f_[pLabels[pi]];
    }
    labelList order(fvert.size());
    sortedOrder(fvert,order);
    scalar f1 = fvert[order.first()];
    scalar f2 = fvert[order.last()];

////    Info << "fvert = " << fvert << ", and order = " << order << endl;
    
    //Handling special case where method is handed an almost full or empty cell
    if (alpha1 < tol)
    {
        return calcSubCell(cellI,f2);
    }
    else if (1-alpha1 < tol)
    {
        return calcSubCell(cellI,f1);
    }
    
    //Finding the two vertices inbetween which the isovalue giving alpha1 lies
    label L1 = 0;
    label L2 = fvert.size()-1;
    scalar a1 = 1;
    scalar a2 = 0;
    scalar L3, f3, a3;
    
    while (L2 - L1 > 1)
    {
        L3 = round(0.5*(L1 + L2));
        f3 = fvert[order[L3]];
        calcSubCell(cellI,f3);
        a3 = VolumeOfFluid();
        if (a3 > alpha1)
        {
            L1 = L3; f1 = f3; a1 = a3; 
        }
        else if (a3 < alpha1)
        {
            L2 = L3; f2 = f3; a2 = a3;
        }        
    }
    
    if (f1 == f2)
    {
        Info << "Warning: f1 = f2." << endl;
        return calcSubCell(cellI,f1);
    }
    
    if (mag(a1-a2) < tol)
    {
        Info << "Warning: mag(a1-a2) < tol for cell " << cellI << endl;
        return calcSubCell(cellI,0.5*(f1 + f2));
    }
    //Now we know that a(f) = alpha1 is to be found on the f interval
    //[f1, f2], i.e. alpha1 will be in the interval [a2,a1]
////    Info << "L1 = " << L1 << ", f1 = " << f1 << ", a1 = " << a1 << endl;
////    Info << "L2 = " << L2 << ", f2 = " << f2  << ", a2 = " << a2 << endl;
    
    
    //Finding coefficients in 3 deg polynomial alpha(f) from 4 solutions
    
    //Finding 2 additional points on 3 deg polynomial
    f3 = f1 + 0.333333333333333333333333333*(f2-f1); 
    calcSubCell(cellI,f3);
    a3 = VolumeOfFluid();

    scalar f4 = f1 + 0.666666666666666666666666667*(f2-f1); 
    calcSubCell(cellI,f4);
    scalar a4 = VolumeOfFluid();
    
    //Building and solving Vandermonde matrix equation
    scalarField a(4), f(4), C(4), fOld(4);
    {
//        f[0] = f1, f[1] = f3, f[2] = f4, f[3] = f2;
        a[0] = a1, a[1] = a3, a[2] = a4, a[3] = a2;
        fOld[0] = f1, fOld[1] = f3, fOld[2] = f4, fOld[3] = f2;
        f[0] = 0, f[1] = (f3-f1)/(f2-f1), f[2] = (f4-f1)/(f2-f1), f[3] = 1;
//        a[0] = 1, a[1] = (a3-a1)/(a1-a2), a[2] = (a4-a1)/(a1-a2), a[3] =1;
        scalarSquareMatrix M(4);
        forAll(f, i)
        {
            forAll(f, j)
            {
                M[i][j] = pow(f[i], 3-j);
            }
        }
        //C holds the 4 polynomial coefficients 
        C = a;
        LUsolve(M, C);
//        Info << "a = " << a << endl;
//        Info << "f = " << f << endl;
//        Info << "C = " << C << endl;
//        Info << "M = " << M << endl;
        if (a[0] < a[1] || a[1] < a[2] || a[2] < a[3])
        {
//            Info << "---------------Cell " << cellI << " with alpha1 = " << alpha1 << endl;
            Info << "Warninig: a is not monotonic: a = " << a << ", f = " << fOld << endl;
/*            
            calcSubCell(cellI,f1);
            a1 = VolumeOfFluid();
            Info << "isoFacePoints for f1 = " << f1 << " and a1 = " << a1 << ": " << endl;
            DynamicList<point> pf = isoFacePoints();
            Info << "p = [" << endl;
            forAll(pf, pi)
            {
                Info << "[" << pf[pi][0] << " " << pf[pi][1] << " " << pf[pi][2] << "];" << endl;
            }
            Info << "];" << endl;
            Info << "isoFaceEdges_ = " << isoFaceEdges_ << endl;

            calcSubCell(cellI,f3);
            a3 = VolumeOfFluid();
            Info << "isoFacePoints for f3 = " << f3 << " and a3 = " << a3 << ": " << endl;
            pf = isoFacePoints();
            Info << "p = [" << endl;
            forAll(pf, pi)
            {
                Info << "[" << pf[pi][0] << " " << pf[pi][1] << " " << pf[pi][2] << "];" << endl;
            }
            Info << "];" << endl;
            Info << "isoFaceEdges_ = " << isoFaceEdges_ << endl;
            
            calcSubCell(cellI,f4);
            a4 = VolumeOfFluid();
            Info << "isoFacePoints for f4 = " << f4 << " and a4 = " << a4 << ": " << endl;
            pf = isoFacePoints();
            Info << "p = [" << endl;
            forAll(pf, pi)
            {
                Info << "[" << pf[pi][0] << " " << pf[pi][1] << " " << pf[pi][2] << "];" << endl;
            }
            Info << "];" << endl;
            Info << "isoFaceEdges_ = " << isoFaceEdges_ << endl;

            calcSubCell(cellI,f2);
            a2 = VolumeOfFluid();
            Info << "isoFacePoints for f2 = " << f2 << " and a2 = " << a2 << ": " << endl;
            pf = isoFacePoints();
            Info << "p = [" << endl;
            forAll(pf, pi)
            {
                Info << "[" << pf[pi][0] << " " << pf[pi][1] << " " << pf[pi][2] << "];" << endl;
            }
            Info << "];" << endl;
            Info << "isoFaceEdges_ = " << isoFaceEdges_ << endl;
            
            Info << "Cell points: " << endl;
            const labelList& pl = mesh_.cellPoints(cellI);
            const pointField& points = mesh_.points();
            Info << "cp = [" << endl;
            forAll(pl, pi)
            {
                Info << "[" << points[pl[pi]][0] << " " << points[pl[pi]][1] << " " << points[pl[pi]][2] << "];" << endl;
            }
            Info << "];" << endl;
            Info << "f = [" << endl;
            forAll(pl, pi)
            {
                Info << f_[pl[pi]] << endl;
            }
            Info << "];" << endl;

            Info << "---------------END-------------" << endl;
        */
        }
    }
    
    //Finding root with Newton method
        
    f3 = f[1]; a3 = a[1];
    label nIter = 0;
    scalar res = mag(a3 - alpha1);
    while (res > tol && nIter < 10*maxIter)
    {
        f3 -= (C[0]*pow3(f3) + C[1]*sqr(f3) + C[2]*f3 + C[3] - alpha1)
            /(3*C[0]*sqr(f3) + 2*C[1]*f3 + C[2]);
        a3 = C[0]*pow3(f3) + C[1]*sqr(f3) + C[2]*f3 + C[3];
        res = mag(a3 - alpha1);
        nIter++;
    }
    //Scaling back to original range
    f3 = f3*(f2 - f1) + f1;

/*    
    if (res > tol)
    {
        Info << "Warning: Leaving Newton method in iter " << nIter 
            << " with f3 = " << f3 << " and a3 = " << a3 << endl;
    }
*/
    //Check result
    calcSubCell(cellI,f3);
    scalar VOF = VolumeOfFluid();
    res = mag(VOF - alpha1);
    
    if (res > tol)
    {
/*
        Info << "Newton obtained f3 = " << f3 << " and a3 = " << a3 
           << " with mag(a3-alpha1) = " << mag(a3-alpha1) 
           << " but calcSubCell(cellI,f3) gives VOF  = " << VOF << endl;
        Info << "M(f)*C = a with " << endl;
        Info << "f_scaled = " << f << endl;
        Info << "f = " << f*(f2 - f1) + f1 << endl;
        Info << "a = " << a << endl;
        Info << "C = " << C << endl;
*/
    }
    else
    {
//        Info << "Newton did the job" << endl;
        return cellStatus_;
    }

    //If tolerance not met use the secant method  with f3 as a hopefully very 
    //good initial guess to crank res the last piece down below tol

    scalar x2 = f3;
    scalar g2 = VOF - alpha1;    
    scalar x1 = max(1e-3*(f2-f1),100*SMALL);
    x1 = max(x1,f1);
    x1 = min(x1,f2);
    calcSubCell(cellI,x1);
    scalar g1 = VolumeOfFluid() - alpha1;

/*
    scalar x1 = f1;
    scalar g1 = a1 - alpha1;    
    scalar x2 = f2;
    scalar g2 = a2 - alpha1;
*/
    nIter = 0;
    scalar g0(0), x0(0);
    while ( res > tol && nIter < maxIter && g1 != g2 )
    {
        x0 = (x2*g1 - x1*g2)/(g1 - g2);
        calcSubCell(cellI,x0);
        g0 = VolumeOfFluid() - alpha1;
        res = mag(g0);
        x2 = x1; g2 = g1;
        x1 = x0; g1 = g0;
        nIter++;
    }
/*
    if (nIter > 0)
    {
        f3 = x0;
        a3 = g0 + alpha1;
    }    
    
    if (res < tol)
    {
        Info << "Bisection finished the job in " << nIter << " iterations." << endl;
    }
    else
    {
        Info << "Warning: Bisection not converged " << endl;
        Info << "Leaving vofCutCell with f3 = " << f3 << " giving a3 = " 
            << a3 << " so alpha1 - a3 = " << alpha1 - a3 << endl;
    }    
*/        
    return cellStatus_;
}




// ************************************************************************* //