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

#ifdef DETAILS2LOG
#define isoDebug(x) x
#else
#define isoDebug(x)
#endif

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
    if ( nPoints == 3 )
    {
        subFaceCentre_ = (1.0/3.0)*(subFacePoints_[0] + subFacePoints_[1]
            + subFacePoints_[2]);
        subFaceArea_ = 0.5*((subFacePoints_[1]
            - subFacePoints_[0])^(subFacePoints_[2] - subFacePoints_[0]));
    }
    else if ( nPoints > 0 )
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
            vector n = (nextPoint
                - subFacePoints_[pi])^(fCentre - subFacePoints_[pi]);
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

    //If vertex values are very close to isoValue lift them slightly to avoid
    //dealing with the many special cases of a face being touched either at a
    //single point, along an edge, or the entire face being on the surface.
    if (mag(f1 - isoValue_) < 10*SMALL)
    {
        f1 += sign(f1 - isoValue_)*10*SMALL;
//        Info << "Warning: adding small number to vertex value of face "
//            << faceI_ << " with owner cell " << mesh_.owner()[faceI_] << endl;
    }

    //Finding cut edges, the point along them where they are cut, and all fully
    //submerged face points.
    forAll(pLabels, pi)
    {
        label pl2 = pLabels[(pi + 1) % nPoints];
        scalar f2 = f[pl2];
        if (mag(f2 - isoValue_) < 10*SMALL)
        {
            f2 += sign(f2 - isoValue_)*10*SMALL;
//            Info << "Warning: adding small number to vertex value of face "
//                << faceI_ << " with owner cell " << mesh_.owner()[faceI_] << endl;
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
                FatalErrorIn
                (
                    "void Foam::isoCutFace::calcSubFace(...)"
                )   << "More than two face cuts for face " << faceI_
                    << abort(FatalError);

                Info << "Warning: More than two face cuts for face " << faceI_ << endl;
                const labelList& fl = mesh_.faces()[faceI_];
                Info << "Face values: f-isoValue = " << endl;
                forAll(fl, fi)
                {
                    Info << f_[fl[fi]]-isoValue_ << " ";
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
        Info << "], " << "firstFulSubPt_: " << firstFullySubmergedPoint_
            << ", nFullySubmergedPoints_: " << nFullySubmergedPoints_
            << " with isoValue_: " << isoValue_ << endl;
*/
    }
    else if (f1 < isoValue_ ) //firstFullySubmergedPoint_ = -1 means no cuttings
    {
        faceStatus_ = 1; //face entirely above isosurface
    }
    //else if (f1 > isoValue_) {face below isosurface, faceStatus_ = -1
    //which is its default value, so no action required here
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

    label pl1 = pLabels[(firstFullySubmergedPoint_
        + nFullySubmergedPoints_ - 1) % nPoints];
    label pl2 = pLabels[(firstFullySubmergedPoint_
        + nFullySubmergedPoints_) % nPoints];

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
    const labelList pLabels = identity(f.size());
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


void Foam::isoCutFace::cutPoints
(
    const pointField& pts,
    const scalarField& f,
    const scalar f0,
    DynamicList<point>& cutPoints
)
{
    isoDebug(Info << "Enter getFaceCutPoints" << endl;)

    const label nPoints = pts.size();
    scalar f1(f[0]);
    forAll(pts,pi)
    {
        label pi2 = (pi+1)%nPoints;
        scalar f2 = f[pi2];
        if ( (f1 < f0 && f2 > f0 ) || (f1 > f0 && f2 < f0) )
        {
            scalar s = (f0-f1)/(f2-f1);
            point pCut = pts[pi] + s*(pts[pi2]-pts[pi]);
            cutPoints.append(pCut);
        }
        else if ( f1 == f0 )
        {
            cutPoints.append(pts[pi]);
        }
        f1 = f2;
    }
}

// ************************************************************************* //
