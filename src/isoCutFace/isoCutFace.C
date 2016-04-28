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
    faceStatus_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::isoCutFace::calcSubFaceCentreAndArea()
{
    
    subFacePoints();

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
}


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
    calcSubFace(f_, pLabels);
    return faceStatus_;
}


Foam::label Foam::isoCutFace::calcSubFace
(
    const scalarField& f,
    const scalar isoValue
)
{    
    clearStorage();
    isoValue_ = isoValue;
    const labelList pLabels = identity(f.size());    
    calcSubFace(f, pLabels);
    return faceStatus_;
}


void Foam::isoCutFace::calcSubFace
(   
    const scalarField& f, 
    const labelList& pLabels
)
{
    const label nPoints = pLabels.size();
    label pl1 = pLabels[0];
    scalar f1 = f[pl1];
    
    if (f1 == isoValue_)
    {
        f1 += 10*SMALL;
    }

    forAll(pLabels, pi)
    {
        label pl2 = pLabels[(pi + 1) % nPoints];
        scalar f2 = f[pl2];
        if (f2 == isoValue_)
        {
            f2 += 10*SMALL;
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
                FatalErrorIn
                (
                    "void Foam::isoCutFace::calcSubFace(...)"
                )   << "More than two face cuts for face " << faceI_
                    << abort(FatalError);
            }
        }
        pl1 = pl2;
        f1 = f2;
    }

    if (firstFullySubmergedPoint_ != -1)
    {
        faceStatus_ = 0;
        Info << "f = [";
        forAll(pLabels, pi)
        {
            scalar f1 = f[pLabels[pi]];
            Info << f1 << ", ";
        }
        Info << "], " << "firstFulSubPt_: " << firstFullySubmergedPoint_ 
            << ", nFullySubmergedPoints_: " << nFullySubmergedPoints_ 
            << " with isoValue_: " << isoValue_ << endl;
    }
    else if (f1 < isoValue_ ) //firstFullySubmergedPoint_ mans no cuttings
    {
        faceStatus_ = 1; //face entirely above isosurface
    } 
    //else if (f1 > isoValue_) {face below isosurface, faceStatus_ = -1
    //which is its default value, so no action required here
}


Foam::point Foam::isoCutFace::subFaceCentre()
{
    calcSubFaceCentreAndArea();
    return subFaceCentre_;
}


Foam::vector Foam::isoCutFace::subFaceArea()
{
    calcSubFaceCentreAndArea();
    return subFaceArea_;
}


Foam::DynamicList<Foam::point> Foam::isoCutFace::subFacePoints()
{
    const labelList& pLabels = mesh_.faces()[faceI_];
    const pointField& points = mesh_.points();    
    subFacePoints(points, pLabels);
    return subFacePoints_;
}


Foam::DynamicList<Foam::point> Foam::isoCutFace::subFacePoints
(
    const pointField& points
)
{        
    const labelList pLabels = identity(points.size());
    subFacePoints(points, pLabels);
    return subFacePoints_;
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


Foam::DynamicList<Foam::point> Foam::isoCutFace::surfacePoints()
{        
    const labelList& pLabels = mesh_.faces()[faceI_];
    const pointField& points = mesh_.points();  
    surfacePoints(points, pLabels);
    return surfacePoints_;
}


Foam::DynamicList<Foam::point> Foam::isoCutFace::surfacePoints
(
    const pointField& points
)
{        
    const labelList pLabels = identity(points.size());
    surfacePoints(points, pLabels);
    return surfacePoints_;
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


void Foam::isoCutFace::clearStorage()
{
    firstEdgeCut_ = -1;
    lastEdgeCut_ = -1;
    firstFullySubmergedPoint_ = -1;
    nFullySubmergedPoints_ = 0;
    faceStatus_ = -1;
    
    subFacePoints_.clear();
    surfacePoints_.clear();
    subFaceCentre_ = vector::zero;
    subFaceArea_ = vector::zero;
}


// ************************************************************************* //