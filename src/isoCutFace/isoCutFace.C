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
    AreaOfFluid_(-1)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::isoCutFace::calcSubFaceCentreAndArea()
{
    
    List<point> p = subFacePoints();

    const label nPoints = p.size();
    // If the face is a triangle, do a direct calculation for efficiency
    // and to avoid round-off error-related problems
    if ( nPoints == 3 )
    {
        subFaceCentre_ = (1.0/3.0)*(p[0] + p[1] + p[2]);
        subFaceArea_ = 0.5*((p[1] - p[0])^(p[2] - p[0]));
    }
    else if ( nPoints > 0 )
    {
        vector sumN = vector::zero;
        scalar sumA = 0.0;
        vector sumAc = vector::zero;

        point fCentre = p[0];
        for (label pi = 1; pi < nPoints; pi++)
        {
            fCentre += p[pi];
        }

        fCentre /= nPoints;

        for (label pi = 0; pi < nPoints; pi++)
        {
            const point& nextPoint = p[(pi + 1) % nPoints];

            vector c = p[pi] + nextPoint + fCentre;
            vector n = (nextPoint - p[pi])^(fCentre - p[pi]);
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
    
    AreaOfFluid_ = mag(subFaceArea_)/mesh_.magSf()[faceI_];
}


Foam::label Foam::isoCutFace::calcSubFace
(
    const label faceI,
    const scalar isoValue
)
{
//    Info << "Entering calcSubFace for face " << faceI << " with isovalue " 
//        << isoValue << endl;
        
    clearStorage();
    faceI_ = faceI;
    isoValue_ = isoValue;
    const labelList& pLabels = mesh_.faces()[faceI_];

    const label nPoints = pLabels.size();
    label pl1 = pLabels[0];
    scalar f1 = f_[pl1];

    forAll(pLabels, pi)
    {
        label pl2 = pLabels[(pi + 1) % nPoints];
        scalar f2 = f_[pl2];

        if (f1 > isoValue_)
        {
            nFullySubmergedPoints_ += 1;
            
            if (f2 <= isoValue_)
            {
                lastEdgeCut_ = (isoValue_ - f1)/(f2 - f1);
            }
        }
        else if (f1 <= isoValue_ && f2 > isoValue_)
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

    if (firstFullySubmergedPoint_!=-1)
    {
        Info << "f = [";
        forAll(pLabels, pi)
        {
            scalar f1 = f_[pLabels[pi]];
            Info << f1 << ", ";
        }
        Info << "], " << "firstFulSubPt_: " << firstFullySubmergedPoint_ 
            << ", nFullySubmergedPoints_: " << nFullySubmergedPoints_ 
            << " with isoValue_: " << isoValue_ << endl;
    }
    return firstFullySubmergedPoint_;
}


Foam::point Foam::isoCutFace::subFaceCentre()
{
    if ( AreaOfFluid_ < 0 )
    {
        calcSubFaceCentreAndArea();
    }

    return subFaceCentre_;
}


Foam::vector Foam::isoCutFace::subFaceArea()
{
    if ( AreaOfFluid_ < 0 )
    {
        calcSubFaceCentreAndArea();
    }

    return subFaceArea_;
}


Foam::List<Foam::point> Foam::isoCutFace::subFacePoints()
{
//    Info << "Enter subFacePoints()" << endl;

    const labelList& pLabels = mesh_.faces()[faceI_];
    const pointField& points = mesh_.points();
    const label nPoints = pLabels.size();
    
    List<point> surfPts = surfacePoints();
    const label nSurfPts = surfPts.size();
    const label nSubFacePoints = nSurfPts + nFullySubmergedPoints_;

    List<point> subPoints(nSubFacePoints);
    forAll(surfPts, pi)
    {
        subPoints[pi] = surfPts[pi];
    }
    
    for(label pi = 0; pi < nFullySubmergedPoints_; pi++)
    {
        subPoints[pi + nSurfPts] = 
            points[pLabels[(firstFullySubmergedPoint_ + pi) % nPoints]];
    }

    return subPoints;
}


Foam::List<Foam::point> Foam::isoCutFace::surfacePoints()
{
//    Info << "Enter surfacePoints()" << endl;
    
    const labelList& pLabels = mesh_.faces()[faceI_];
    const pointField& points = mesh_.points();
    const label nPoints = pLabels.size();
    
    label pl1 = pLabels[(firstFullySubmergedPoint_ 
        + nFullySubmergedPoints_ - 1) % nPoints];
    label pl2 = pLabels[(firstFullySubmergedPoint_ 
        + nFullySubmergedPoints_) % nPoints];

    List<point> surfacePoints(2);
    surfacePoints[0] = points[pl1] + lastEdgeCut_*(points[pl2] - points[pl1]);
    
    pl1 = pLabels[(firstFullySubmergedPoint_ - 1 + nPoints) % nPoints];
    pl2 = pLabels[firstFullySubmergedPoint_];

    surfacePoints[1] = points[pl1] + firstEdgeCut_*(points[pl2] - points[pl1]);
//    Info << "surfPoints for face" << faceI_ << ": " << surfacePoints << endl;

    return surfacePoints;
}


Foam::scalar Foam::isoCutFace::areaOfFluid()
{
    if ( AreaOfFluid_ < 0 )
    {
        calcSubFaceCentreAndArea();
    }
    
    return AreaOfFluid_;
}


void Foam::isoCutFace::clearStorage()
{
    firstEdgeCut_ = -1;
    lastEdgeCut_ = -1;
    firstFullySubmergedPoint_ = -1;
    nFullySubmergedPoints_ = 0;

    subFaceCentre_ = vector::zero;
    subFaceArea_ = vector::zero;
    AreaOfFluid_ = -1;
}


// ************************************************************************* //