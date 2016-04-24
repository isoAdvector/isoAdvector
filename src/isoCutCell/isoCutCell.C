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
    cellStatus_(-1)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::isoCutCell::calcSubCellCentreAndVolume()
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
}


void Foam::isoCutCell::calcIsoFaceCentreAndNormal()
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
    fCentre /= nEdgePoints;
    
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
}


void Foam::isoCutCell::calcIsoFacePointsFromEdges()
{
    
    forAll(isoFaceEdges_, ei)
    {

        const DynamicList<point>& edgePoints = isoFaceEdges_[ei];
        const label nPoints = edgePoints.size();
        label edgeOrientation = 0;
        if (nPoints == 2)
        {
            const point& P1 = edgePoints[0];
            const point& P2 = edgePoints[1];
            edgeOrientation = 
                sign(isoFaceArea_ & ((P1-isoFaceCentre_)^(P2-isoFaceCentre_)));
            if ( edgeOrientation == -1 )
            {
                isoFacePoints_.append(P1);
            }
            else
            {
                isoFacePoints_.append(P2);                
            }
        }
        else if (nPoints == 1)
        {
            isoFacePoints_.append(edgePoints[0]);
        }
        else if (nPoints > 2)
        {
            Info << "Warning: A face of cell " << cellI_ << " has " << nPoints 
                << " cut points." << endl;
        }
    }
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
    const labelList& faces = mesh_.faces()[cellI];
    
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
        calcIsoFaceCentreAndNormal();
    }
    else if (fullySubFaces_.size() == 0) //cell fully above isosurface
    {
        cellStatus_ = 1;
    }
    //else all faces are below isosurface and cellStatus_ = -1,
    //which is its default value, so no action required here

    return cellStatus_;
}


Foam::point Foam::isoCutCell::subCellCentre()
{
    if ( VOF_ < -9 )
    {
        calcSubCellCentreAndVolume();
    }

    return subCellCentre_;
}


Foam::scalar Foam::isoCutCell::subCellVolume()
{
    if ( VOF_ < -9 )
    {
        calcSubCellCentreAndVolume();
    }

    return subCellVolume_;
}


Foam::DynamicList<Foam::point> Foam::isoCutCell::isoFacePoints()
{
    return isoFacePoints_;
}


Foam::point Foam::isoCutCell::isoFaceCentre()
{
    return isoFaceCentre_;
}

Foam::vector Foam::isoCutCell::isoFaceArea()
{    
    return isoFaceArea_;
}
        

Foam::scalar Foam::isoCutCell::VolumeOfFluid()
{
    if ( VOF_ < -9 )
    {
        calcSubCellCentreAndVolume();
    }
    
    return VOF_;
}


void Foam::isoCutCell::clearStorage()
{
    isoCutFaces_.clear();
    isoCutFacePoints_.clear();
    isoCutFaceCentres_.clear();
    isoCutFaceAreas_.clear();
    isoFacePoints_.clear();
    isoFaceEdges_.celar();
    fullySubFaces_.clear();
    VOF_ = -10;
    cellStatus_ = -1;
}


// ************************************************************************* //