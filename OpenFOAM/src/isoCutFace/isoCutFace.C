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

#include "isoCutFace.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isoCutFace::isoCutFace
(
    const fvMesh& mesh,
    scalarField& f
)
:
    mesh_(mesh),
    faceI_(-1),
    f_(f),
    isoValue_(0),
    firstEdgeCut_(-1),
    lastEdgeCut_(-1),
    firstFullySubmergedPoint_(-1),
    nFullySubmergedPoints_(0),
    subFaceCentre_(vector::zero),
    subFaceArea_(vector::zero),
    subFacePoints_(10),
    surfacePoints_(4),
    faceStatus_(-1),
    subFaceCentreAndAreaIsCalculated_(false)
{
    clearStorage();
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::isoCutFace::calcSubFaceCentreAndArea()
{
    const label nPoints = subFacePoints_.size();

    // If the face is a triangle, do a direct calculation for efficiency
    // and to avoid round-off error-related problems
    if (nPoints == 3)
    {
        subFaceCentre_ = (1.0/3.0)*(subFacePoints_[0] + subFacePoints_[1]
            + subFacePoints_[2]);
        subFaceArea_ = 0.5*((subFacePoints_[1]
            - subFacePoints_[0])^(subFacePoints_[2] - subFacePoints_[0]));
    }
    else if (nPoints > 0)
    {
        vector sumN = vector::zero;
        scalar sumA = 0.0;
        vector sumAc = vector::zero;

        point fCentre = subFacePoints_[0];
        for (label pi = 1; pi < nPoints; pi++)
        {
            fCentre += subFacePoints_[pi];
        }

        fCentre /= nPoints;

        for (label pi = 0; pi < nPoints; pi++)
        {
            const point& nextPoint = subFacePoints_[(pi + 1) % nPoints];

            vector c = subFacePoints_[pi] + nextPoint + fCentre;
            vector n =
                (nextPoint - subFacePoints_[pi])^(fCentre - subFacePoints_[pi]);
            scalar a = mag(n);

            sumN += n;
            sumA += a;
            sumAc += a*c;
        }

        // This is to deal with zero-area faces. Mark very small faces
        // to be detected in e.g., processorPolyPatch.
        if (sumA < ROOTVSMALL)
        {
            subFaceCentre_ = fCentre;
            subFaceArea_ = vector::zero;
        }
        else
        {
            subFaceCentre_ = (1.0/3.0)*sumAc/sumA;
            subFaceArea_ = 0.5*sumN;
        }
    }

    subFaceCentreAndAreaIsCalculated_ = true;
}


void Foam::isoCutFace::calcSubFace
(
    const pointField& points,
    const scalarField& f,
    const labelList& pLabels
)
{
    const label nPoints = pLabels.size();
    label pl1 = pLabels[0];
    scalar f1 = f[pl1];

    // If vertex values are very close to isoValue lift them slightly to avoid
    // dealing with the many special cases of a face being touched either at a
    // single point, along an edge, or the entire face being on the surface.
    if (mag(f1 - isoValue_) < 10*SMALL)
    {
        f1 += sign(f1 - isoValue_)*10*SMALL;
//        Info<< "Warning: adding small number to vertex value of face "
//            << faceI_ << " with owner cell " << mesh_.faceOwner()[faceI_]
//            << endl;
    }

    // Finding cut edges, the point along them where they are cut, and all fully
    // submerged face points.
    forAll(pLabels, pi)
    {
        label pl2 = pLabels[(pi + 1) % nPoints];
        scalar f2 = f[pl2];
        if (mag(f2 - isoValue_) < 10*SMALL)
        {
            f2 += sign(f2 - isoValue_)*10*SMALL;
//            Info<< "Warning: adding small number to vertex value of face "
//                << faceI_ << " with owner cell " << mesh_.faceOwner()[faceI_]
//                << endl;
        }

        if (f1 > isoValue_)
        {
            nFullySubmergedPoints_ += 1;

            if (f2 < isoValue_)
            {
                lastEdgeCut_ = (isoValue_ - f1)/(f2 - f1);
            }
        }
        else if (f1 < isoValue_ && f2 > isoValue_)
        {
            if (firstFullySubmergedPoint_ == -1)
            {
                firstFullySubmergedPoint_ = (pi + 1) % nPoints;

                firstEdgeCut_ = (isoValue_ - f1)/(f2 - f1);
            }
            else
            {
/*
                FatalErrorInFunction
                    << "More than two face cuts for face " << faceI_
                    << abort(FatalError);

                Info<< "Warning: More than two face cuts for face " << faceI_
                    << endl;
                const labelList& fl = mesh_.faces()[faceI_];
                Info << "Face values: f-isoValue = " << endl;
                forAll(fl, fi)
                {
                    Info << f_[fl[fi]] - isoValue_ << " ";
                }
                Info << " " << endl;
*/
            }
        }
        pl1 = pl2;
        f1 = f2;
    }

    if (firstFullySubmergedPoint_ != -1)
    {
        faceStatus_ = 0;
        subFacePoints(points, pLabels);
/*
        Info << "f = [";
        forAll(pLabels, pi)
        {
            scalar f1 = f[pLabels[pi]];
            Info << f1 << ", ";
        }
        Info<< "], " << "firstFulSubPt_: " << firstFullySubmergedPoint_
            << ", nFullySubmergedPoints_: " << nFullySubmergedPoints_
            << " with isoValue_: " << isoValue_ << endl;
*/
    }
    else if (f1 < isoValue_) // firstFullySubmergedPoint_ = -1 means no cuttings
    {
        faceStatus_ = 1; // face entirely above isosurface
    }
    // else if (f1 > isoValue_) {face below isosurface, faceStatus_ = -1
    // which is its default value, so no action required here
}


void Foam::isoCutFace::subFacePoints
(
    const pointField& points,
    const labelList& pLabels
)
{
    const label nPoints = pLabels.size();

    surfacePoints(points, pLabels);

    forAll(surfacePoints_, pi)
    {
        subFacePoints_.append(surfacePoints_[pi]);
    }

    for(label pi = 0; pi < nFullySubmergedPoints_; pi++)
    {
        subFacePoints_.append
        (
            points[pLabels[(firstFullySubmergedPoint_ + pi) % nPoints]]
        );
    }
}


void Foam::isoCutFace::surfacePoints
(
    const pointField& points,
    const labelList& pLabels
)
{
    const label nPoints = pLabels.size();

    label pl1 =
        pLabels[(firstFullySubmergedPoint_ + nFullySubmergedPoints_ - 1) % nPoints];

    label pl2 =
        pLabels[(firstFullySubmergedPoint_ + nFullySubmergedPoints_) % nPoints];

    surfacePoints_.append
    (
        points[pl1] + lastEdgeCut_*(points[pl2] - points[pl1])
    );

    pl1 = pLabels[(firstFullySubmergedPoint_ - 1 + nPoints) % nPoints];
    pl2 = pLabels[firstFullySubmergedPoint_];

    surfacePoints_.append
    (
        points[pl1] + firstEdgeCut_*(points[pl2] - points[pl1])
    );
}


// * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //


Foam::label Foam::isoCutFace::calcSubFace
(
    const label faceI,
    const scalar isoValue
)
{
    clearStorage();
    faceI_ = faceI;
    isoValue_ = isoValue;
    const labelList& pLabels = mesh_.faces()[faceI_];
    const pointField& points = mesh_.points();
    calcSubFace(points, f_, pLabels);
    return faceStatus_;
}


Foam::label Foam::isoCutFace::calcSubFace
(
    const pointField& points,
    const scalarField& f,
    const scalar isoValue
)
{
    clearStorage();
    isoValue_ = isoValue;
    const labelList pLabels(identity(f.size()));
    calcSubFace(points, f, pLabels);
    return faceStatus_;
}


Foam::point Foam::isoCutFace::subFaceCentre()
{
    if (!subFaceCentreAndAreaIsCalculated_)
    {
        calcSubFaceCentreAndArea();
    }
    return subFaceCentre_;
}


Foam::vector Foam::isoCutFace::subFaceArea()
{
    if (!subFaceCentreAndAreaIsCalculated_)
    {
        calcSubFaceCentreAndArea();
    }
    return subFaceArea_;
}


Foam::DynamicList<Foam::point> Foam::isoCutFace::subFacePoints()
{
    return subFacePoints_;
}


Foam::DynamicList<Foam::point> Foam::isoCutFace::surfacePoints()
{
    return surfacePoints_;
}


void Foam::isoCutFace::clearStorage()
{
    faceI_ = -1;
    isoValue_ = 0;
    firstEdgeCut_ = -1;
    lastEdgeCut_ = -1;
    firstFullySubmergedPoint_ = -1;
    nFullySubmergedPoints_ = 0;
    subFaceCentre_ = vector::zero;
    subFaceArea_ = vector::zero;
    subFacePoints_.clear();
    surfacePoints_.clear();
    faceStatus_ = -1;
    subFaceCentreAndAreaIsCalculated_ = false;
}


Foam::scalar Foam::isoCutFace::timeIntegratedArea
(
    const pointField& fPts,
    const scalarField& pTimes,
    const scalar dt,
    const scalar magSf,
    const scalar Un0
)
{

    // Initialise time integrated area returned by this function
    scalar tIntArea = 0.0;

    // Finding ordering of vertex points
    labelList order(pTimes.size());
    sortedOrder(pTimes, order);
    const scalar firstTime = pTimes[order.first()];
    const scalar lastTime = pTimes[order.last()];
    
    // Times smaller than tSmall are regarded as 0
    const scalar tSmall = 1e-8*max(lastTime - firstTime, dt);

    // Making sorted list of the labels of those vertex times that are within 
    // the integration interval [0, dt]
    DynamicList<label> cutVertexLabels(pTimes.size());
//    scalar prevTime = 0;
    forAll(order, ti)
    {
        const scalar newTime = pTimes[order[ti]];
//        if (newTime >  prevTime + tSmall &&  newTime < dt - tSmall)
        if (newTime >  0 &&  newTime < dt)
        {
            cutVertexLabels.append(order[ti]);
        }
//        prevTime = newTime;
    }
    
    // Dealing with case where face is not cut by surface during time interval
    // [0,dt] because face was already passed by surface
    if (lastTime < tSmall)
    {
        // If all face cuttings were in the past and cell is filling up (Un0>0)
        // then face must be full during whole time interval
        tIntArea = magSf*dt*pos(Un0);
        return tIntArea;
    }

    // Dealing with case where face is not cut by surface during time interval
    // [0, dt] because dt is too small for surface to reach closest face point
    if (firstTime > dt - tSmall)
    {
        // If all cuttings are in the future but non of them within [0,dt] then
        // if cell is filling up (Un0 > 0) face must be empty during whole time
        // interval
        tIntArea = magSf*dt*(1 - pos(Un0));
        return tIntArea;
    }

    // If we reach this point in the code some part of the face will be swept
    // during [tSmall, dt-tSmall]. However, it may be the case that there are no
    // vertex times within the interval. This will happen sometimes for small 
    // time steps where both the initial and the final face-interface 
    // intersection line (FIIL) will be along the same two edges.
    
    // Face-interface intersection line (FIIL) to be swept across face
    DynamicList<point> FIIL(3);
    // Counter for traversing cut vertices
    label nextVertexLabel = 0;
    // First time in sub time intervals
    scalar tOld = 0;
    // Submerged area at beginning of each sub time interval time
    scalar initialArea = 0.0;
    
    // Special treatment of first sub time interval
    if (firstTime > 0)
    {
        // If firstTime > 0 the face is uncut in the time interval 
        // [0, firstTime] and hence fully submerged in fluid A or B. 
        // If Un0 > 0 cell is filling up and it must initially be empty.
        // If Un0 < 0 cell must initially be full(y immersed in fluid A).
        tOld = firstTime;
        initialArea = magSf*(1.0 - pos(Un0));
        tIntArea = initialArea*tOld;
        cutPoints(fPts, pTimes, pTimes[cutVertexLabels[0]], FIIL);
        nextVertexLabel++;
    }
    else
    {
        // If firstTime <= 0 then face is initially cut and we must
        // calculate the initial submerged area and FIIL:
        calcSubFace(fPts, -sign(Un0)*pTimes, 0.0);
        initialArea = mag(subFaceArea());
        cutPoints(fPts, pTimes, 0.0, FIIL);
    }

    // Now we calculate and add contributions to the time integrated area from 
    // quadrilaterals spanned by consecutive FIIL's up to the last vertex hit
    // in the time interval [0, dt].

    while (nextVertexLabel < cutVertexLabels.size()-1)
    {
        const label newLabel = cutVertexLabels[nextVertexLabel];
        const scalar tNew = pTimes[newLabel];
        if (tNew - tOld > tSmall)
        {
            // New face-interface intersection line
            DynamicList<point> newFIIL(3);
            cutPoints(fPts, pTimes, tNew, newFIIL);
        
            // quadrilateral area coefficients
            scalar alpha = 0, beta = 0;
            quadAreaCoeffs(FIIL, newFIIL, alpha, beta);
            // Integration of area(t) = A*t^2+B*t from t = 0 to 1
            tIntArea += (tNew - tOld)*
                (initialArea + sign(Un0)*(alpha/3.0 + 0.5*beta));
            // Adding quad area to submerged area
            initialArea += sign(Un0)*(alpha + beta);
            FIIL = newFIIL;
        }
        
        tOld = tNew;
        nextVertexLabel++;
    }
    
    // Now, if dt > lastTime the FIIL will leave the face or else it 
    // stopped at the last vertex hit within the time interval [0, dt].
    // In the former case we must add a contribution to tIntArea from the last
    // sub time interval, [lastTime, dt], where pure fluid A or B is 
    // fluxed through the face. In the latter case we must add the final 
    // contribution, where the FIIL sweeps the area from the last hit vertex to
    // its position on the face at time dt.
    
    if (dt > lastTime)
    {
        tIntArea += magSf*(dt - lastTime)*pos(Un0);
    }
    else
    {
        // New face-interface intersection line
        DynamicList<point> newFIIL(3);
        cutPoints(fPts, pTimes, dt, newFIIL);        
        // quadrilateral area coefficients
        scalar alpha = 0, beta = 0;
        quadAreaCoeffs(FIIL, newFIIL, alpha, beta);
        // Integration of area(t) = A*t^2+B*t from t = 0 to 1
        const scalar integratedQuadArea = sign(Un0)*(alpha/3.0 + 0.5*beta);
        tIntArea += (dt - tOld)*(initialArea + integratedQuadArea);
    }
    
    return tIntArea;
}

/*
void Foam::isoCutFace::cutPoints
(
    const pointField& pts,
    const scalarField& f,
    const scalar f0,
    DynamicList<point>& cutPoints
)
{
    const label nPoints = pts.size();
    scalar f1(f[0]);
    forAll(pts, pi)
    {
        label pi2 = (pi + 1) % nPoints;
        scalar f2 = f[pi2];
        if ((f1 < f0 && f2 > f0) || (f1 > f0 && f2 < f0))
        {
            const scalar s = (f0 - f1)/(f2 - f1);
            cutPoints.append(pts[pi] + s*(pts[pi2] - pts[pi]));
        }
        else if (f1 == f0)
        {
            cutPoints.append(pts[pi]);
        }
        f1 = f2;
    }
}
*/

void Foam::isoCutFace::cutPoints
(
    const pointField& pts,
    const scalarField& f,
    const scalar f0,
    DynamicList<point>& cutPoints
)
{
    const label nPoints = pts.size();
    
    const scalar eps = 10*SMALL;
    forAll(f, L1)
    {
        const label L2 = (L1 + 1) % nPoints;
        if (f[L1] < f0 && f[L2] > f0)
        {
            if (f[L2] - f[L1] > eps)
            {
                const scalar s = (f0 - f[L1])/(f[L2] - f[L1]);
                cutPoints.append(pts[L1] + s*(pts[L2] - pts[L1]));
            }
            else
            {
                cutPoints.append(pts[L1]);
            }
        }
        else if (f[L1] > f0 && f[L2] < f0)
        {
            if (f[L1] - f[L2] > eps)
            {
                const scalar s = (f0 - f[L1])/(f[L2] - f[L1]);
                cutPoints.append(pts[L1] + s*(pts[L2] - pts[L1]));
            }
            else
            {
                cutPoints.append(pts[L1]);
            }
        }
        else if (f[L1] == f0)
        {
            cutPoints.append(pts[L1]);
        }
    }
}


void Foam::isoCutFace::quadAreaCoeffs
(
    const DynamicList<point>& pf0,
    const DynamicList<point>& pf1,
    scalar& alpha,
    scalar& beta
) const
{

    const label np0 = pf0.size();
    const label np1 = pf1.size();

    alpha = 0.0;
    beta = 0.0;

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
            if (np0 != 1)
            {
                WarningInFunction
                    << "Vertex face was cut at pf0 = " << pf0 << endl;
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
                WarningInFunction
                    << "Vertex face was cut at pf1 = " << pf1 << endl;
            }
        }

        // Defining local coordinates for area integral calculation
        vector xhat = B - A;
        xhat /= mag(xhat);
        vector yhat = D - A;
        yhat -= ((yhat & xhat)*xhat);
        yhat /= mag(yhat);

//            Info<< "xhat = " << xhat << ", yhat = " << yhat << ", zhat = "
//                << zhat << ". x.x = " << (xhat & xhat) << ", y.y = "
//                << (yhat & yhat) <<", z.z = " << (zhat & zhat) << ", x.y = "
//                << (xhat & yhat) << ", x.z = " << (xhat & zhat) << ", y.z = "
//                << (yhat & zhat) << endl;

        // Swapping pf1 points if pf0 and pf1 point in same general direction
        // (because we want a quadrilateral ABCD where pf0 = AB and pf1 = CD)
        if (((B - A) & (D - C)) > 0)
        {
            vector tmp = D;
            D = C;
            C = tmp;
        }

        const scalar Bx = mag(B - A);
        const scalar Cx = (C - A) & xhat;
        const scalar Cy = mag((C - A) & yhat);
        const scalar Dx = (D - A) & xhat;
        const scalar Dy = mag((D - A) & yhat);

//      area = ((Cx-Bx)*Dy-Dx*Cy)/6.0 + 0.25*Bx*(Dy+Cy);
        alpha = 0.5*((Cx-Bx)*Dy-Dx*Cy);
        beta = 0.5*Bx*(Dy+Cy);
//        quadArea = alpha + beta;

//        intQuadArea = alpha/3.0 + 0.5*beta;

//         Info<< "Bx = " << Bx << ", Cx = " << Cx << ", Cy = " << Cy
//             << ", Dx = " << Dx << ", Dy = " << Dy << ", alpha = " << alpha
//             << ", beta = " << beta << endl;

        // area(t) = A*t^2 + B*t
        // integratedArea = A/3 + B/2
    }
    else
    {
        Info<< "Vertex face was cut at " << pf0
            << " and at " << pf1 << " by " << endl;
    }
}


// ************************************************************************* //
