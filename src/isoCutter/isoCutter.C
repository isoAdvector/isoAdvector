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

#include "isoCutter.H"
#include "volPointInterpolation.H" //Is this used at all?
#include "interpolationCellPoint.H"
//#include "isoSubCell.H"

#ifdef ISODEBUG
#define isoDebug(x) x
#else
#define isoDebug(x) 
#endif 

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isoCutter::isoCutter
(
    const fvMesh& mesh,
    const scalarField& f//,
//    const scalar& f0
)
:
    mesh_(mesh),
    f_(f)//,
//    f0_(f0)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::isoCutter::vofCutCell
(
    const label ci,
    const scalar alpha1,
    const scalar tol,
    const label maxIter,
    scalar& f0,
    vector& subCellCtr
)
{
////    Info << "Enter vofCutCell for cell " << ci << " which has alpha1 = " << alpha1 << endl;
    scalar fMin(GREAT), fMax(-GREAT);
    const labelList& pLabels = mesh_.cellPoints(ci);
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
//    Info << "fMin = " << fMin << " fMax = " << fMax << endl;

    //Initial guess of isovalue
    scalar aMin(0), aMax(1), alpha0;
    f0 = (alpha1 - aMin)/(aMax-aMin)*(fMin-fMax) + fMax;
//  f0 = 0.5*(fMin + fMax);
    subCellFraction(ci, f0, alpha0, subCellCtr);
//    Info << "f0 = " << f0 << " gives alpha = " << alpha0 << endl;

    //Bisection method to find root
    //Since function is monotonically increasing and only piecewise smooth derivative based root finding algorithms should not be used here.
    label nIter = 0;
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
//      f0 = (alpha1 - aMin)/(aMax-aMin)*(fMin-fMax) + fMax; //This does an extremely poor job in narrowing in the interval - especially for almost filled cells
        subCellFraction(ci, f0, alpha0, subCellCtr);
//        Info << nIter << ": f0 = " << f0 << " gives alpha = " << alpha0 << endl;
        nIter++;
    }
//    Info << nIter-1 << ": f0 = " << f0 << " gives alpha = " << alpha0 << endl;
}

void Foam::isoCutter::isoCutCell
(
    const label cellI,
    const scalar f0,
    DynamicList<label>& cutFaces,
    DynamicList<label>& cutEdges,
    DynamicList<scalar>& cutPoints,
    bool& fullySubmerged
)
{
    const cellList& cells = mesh_.cells(); //It is important for the efficiency of the code to verify that this call does not generate the cell faces for all cells each time it is called but only the first time!
    labelList cellFaces = cells[cellI];

    label face, edge;
    scalar cutPoint;
    findACutFaceEdgePair(cellI,f0,face,edge,cutPoint,fullySubmerged);
    if ( edge != -1 ) //Only entered if a cut face is found for this cell in which case we should be able to go all the way around
    {
        label firstFace = face;
        while ( true )
        {
//            Info << "Face: " << cellFaceLabels[face] << ", edge: " << edge << ", cutPoint: " << cutPoint << endl;
            cutFaces.append(cellFaces[face]);
            cutEdges.append(edge);
            cutPoints.append(cutPoint);
            findNextCut(cellI, f0, face, edge, cutPoint);
            if ( face == firstFace )
            {
                break;
            }
        }
        cutFaces.shrink();
        cutEdges.shrink();
        cutPoints.shrink();
    }
}


void Foam::isoCutter::findACutFaceEdgePair
(
    const label cellI,
    const scalar f0,
    label& cutFace,
    label& cutEdge,
    scalar& cutPoint,
    bool& fullySubmerged
)
{
    cutEdge = -1;
    cutFace = -1;
    cutPoint = -1;
    fullySubmerged = true;

    const faceList& faces = mesh_.faces();
    const labelList& cellFaces = mesh_.cells()[cellI];

    labelList pLabels;
    forAll(cellFaces,fi)
    {
        pLabels = faces[cellFaces[fi]];
        label nPoints = pLabels.size();
        scalar f1 = f_[pLabels[0]];
        forAll(pLabels,pi)
        {
            scalar f2 = f_[pLabels[(pi+1) % nPoints]];
            if ( min(f1,f2) <= f0 )
            {
                fullySubmerged = false;

                if ( max(f1,f2) > f0 )
                {
                    cutEdge = pi;
                    cutFace = fi;
                    cutPoint = (f0-f1)/(f2-f1);
                    break;
                }
            }
            f1 = f2;
        }
        if ( cutEdge != -1 )
        {
            break;
        }
    }
}


void Foam::isoCutter::findNextCut
(
    const label cellI,
    const scalar f0,
    label& cellFace,
    label& edge,
    scalar& cutPoint
)
{

    const label oldEdge = edge;
    cutPoint = -1;
    const faceList& faces = mesh_.faces();
    const labelList& cellFaces = mesh_.cells()[cellI];

    labelList pLabels = faces[cellFaces[cellFace]];
    label nPoints = pLabels.size();
    edge = (edge + 1) % nPoints;
    scalar f1 = f_[pLabels[edge]];
    while ( edge != oldEdge )
    {
        scalar f2 = f_[pLabels[(edge + 1) % nPoints]];
        if ( min(f1,f2) <= f0 && max(f1,f2) > f0 ) //By convention a face is partially submerged if at least one point is above or on the surface and at least one point is strictly below the surface
        {
            cutPoint = (f0-f1)/(f2-f1);
            label p1 = pLabels[edge];
            otherEdgeFace(cellI, cellFace, edge);
            pLabels = faces[cellFaces[cellFace]];
            if (  pLabels[edge] != p1 )
            {
                cutPoint = 1-cutPoint;
            }
            break;
        } else
        {
            edge = (edge + 1) % nPoints;
            f1 = f2;
        }
    }
}


void Foam::isoCutter::otherEdgeFace
(
    const label cellI,
    label& cellFace,
    label& edge
)
{
    const faceList& faces = mesh_.faces();
    const label oldFace = cellFace;
    const labelList& cellFaces = mesh_.cells()[cellI];

    labelList pLabels = faces[cellFaces[oldFace]];
    label p1 = pLabels[edge];
    label nPoints = pLabels.size();
    label p2 = pLabels[(edge+1) % nPoints];

    cellFace = -1;
    edge = -1;

    forAll(cellFaces,fi)
    {
        if ( fi != oldFace )
        {
            pLabels = faces[cellFaces[fi]];
            nPoints = pLabels.size();
            forAll(pLabels,pi)
            {
                label pMin = min(pLabels[pi],pLabels[(pi+1) % nPoints]);
                label pMax = max(pLabels[pi],pLabels[(pi+1) % nPoints]);
                if ( pMin == min(p1,p2) && pMax == max(p1,p2) )
                {
                    edge = pi;
                    break;
                }
            }
            if ( edge != -1 )
            {
                cellFace = fi;
                break;
            }
        }
    }
}


void Foam::isoCutter::getIsoFace
(
    const DynamicList<label>& cutFaces,
    const DynamicList<label>& cutEdges,
    const DynamicList<scalar>& cutPoints,
    DynamicList<point>& isoPoints
)
{
    const faceList& faces = mesh_.faces();
    const pointField& points = mesh_.points();
    bool findDuplicates = false;

    DynamicList<point> p(cutFaces.size());
    forAll(cutFaces,fi)
    {
        labelList pLabels = faces[cutFaces[fi]];
        label l0 = pLabels[cutEdges[fi]];
        label nPoints = pLabels.size();
        label l1 = pLabels[(cutEdges[fi]+1) % nPoints];
        p.append(points[l0] + cutPoints[fi]*(points[l1]-points[l0]));
        if ( cutPoints[fi] == 0 || cutPoints[fi] == 1 )
        {
//            Info << "s = " << cutPoints[fi] << endl;
            findDuplicates = true;
        }
    }

    //Removing duplicates
    if ( findDuplicates )
    {
        DynamicList<point> pCleaned(p.size());
        label nPoints = p.size();
        forAll(p,pi)
        {
            if ( mag(p[pi]-p[(pi+1)%nPoints]) > VSMALL )
            {
                pCleaned.append(p[pi]);
            }
        }
        if ( pCleaned.size() < 1 && p.size() > 0 )
        {
            pCleaned.append(p[0]);
        }
        isoPoints = pCleaned;
    }
    else
    {
        isoPoints = p;
    }
}

bool Foam::isoCutter::getSubFace
(
    const label faceLabel,
    const scalar f0,
    DynamicList<point>& partSubFacePts
)
{
    const faceList& faces = mesh_.faces();
    const pointField& points = mesh_.points();
    bool fullySubmerged = true;

    const labelList pLabels = faces[faceLabel];
    const label nPoints = pLabels.size();
    label pl1 = pLabels[0];
    forAll(pLabels,pi)
    {
        label pl2 = pLabels[(pi+1)%nPoints];
        scalar f1(f_[pl1]), f2(f_[pl2]);
        if (f1 >= f0)
        {
            partSubFacePts.append(points[pl1]);
            if ( f2 < f0 && f1 > f0 )
            {
                scalar s = (f0-f1)/(f2-f1);
//                Info << "Face " << faceLabel << ", edge " << pi << ": s = " << s << endl;
                point pCut = points[pl1] + s*(points[pl2]-points[pl1]);
                partSubFacePts.append(pCut);
            }
        }
        else if (f1 < f0)
        {
            fullySubmerged = false;
            if (f2 > f0)
            {
                scalar s = (f0-f1)/(f2-f1);
//                Info << "Face " << faceLabel << ", edge " << pi << ": s = " << s << endl;
                point pCut = points[pl1] + s*(points[pl2]-points[pl1]);
                partSubFacePts.append(pCut);
            }
        }
        pl1 = pl2;
    }

    return fullySubmerged;
}


Foam::scalar Foam::isoCutter::getSubFaceFraction
(
    const label faceLabel,
    const scalarField& f,
    const scalar f0
)
{
    const faceList& faces = mesh_.faces();
    const pointField& points = mesh_.points();
    bool fullySubmerged = true;
    const labelList pLabels = faces[faceLabel];
    const label nPoints = pLabels.size();

    DynamicList<point> partSubFacePts(nPoints);

    label pl1 = pLabels[0];
    forAll(pLabels,pi)
    {
        label pl2 = pLabels[(pi+1)%nPoints];
        scalar f1(f[pi]), f2(f[(pi+1)%nPoints]);
        if (f1 >= f0)
        {
            partSubFacePts.append(points[pl1]);
            if ( f2 < f0 && f1 > f0 )
            {
                scalar s = (f0-f1)/(f2-f1);
//                Info << "Face " << faceLabel << ", edge " << pi << ": s = " << s << endl;
                point pCut = points[pl1] + s*(points[pl2]-points[pl1]);
                partSubFacePts.append(pCut);
            }
        }
        else if (f1 < f0)
        {
            fullySubmerged = false;
            if (f2 > f0)
            {
                scalar s = (f0-f1)/(f2-f1);
//                Info << "Face " << faceLabel << ", edge " << pi << ": s = " << s << endl;
                point pCut = points[pl1] + s*(points[pl2]-points[pl1]);
                partSubFacePts.append(pCut);
            }
        }
        pl1 = pl2;
    }
    partSubFacePts.shrink();

    //Calculating subface fraction
    scalar alphaf = 0.0;

    if (fullySubmerged)
    {
        alphaf = 1;
    }
    else if (partSubFacePts.size() != 0)
    {
        vector fCtr, fArea;
        makeFaceCentreAndArea(partSubFacePts, fCtr, fArea);
        alphaf = mag(fArea)/mag(mesh_.Sf()[faceLabel]);
    }

    return alphaf;
}


void Foam::isoCutter::getFaceCutPoints
(
    const label fLabel,
    const scalar f0,
    DynamicList<point>& cutPoints
)
{
    const faceList& faces = mesh_.faces();
    const pointField& points = mesh_.points();

    const labelList pLabels = faces[fLabel];
    const label nPoints = pLabels.size();
    label pl1 = pLabels[0];
    scalar f1(f_[pl1]);
//    if (mag(f1-f0) < 1e-10)
//    {
//        f1 = f0;
//    }
    Info << "face values differences from f0: ";
    forAll(pLabels,pi)
    {
        label pl2 = pLabels[(pi+1)%nPoints];
        scalar f2(f_[pl2]);
        if (mag(f2-f0) < 1e-12)
        {
            f2 = f0;
        }
        Info << " " << f1-f0;
        bool edgeIsCut = min(f1,f2) < f0 && max(f1,f2) > f0;
        if ( edgeIsCut )
        {
            scalar s = (f0-f1)/(f2-f1);
            point pCut = points[pl1] + s*(points[pl2]-points[pl1]);
            cutPoints.append(pCut);
        }
        else if ( f1 == f0 )
        {
            cutPoints.append(points[pl1]);
        }
        pl1 = pl2;
        f1 = f2;
    }
    Info << "." << endl;
}

void Foam::isoCutter::getFaceCutPoints
(
    const label fLabel,
    const scalarField& f,
    const scalar f0,
    DynamicList<point>& cutPoints
)
{
	isoDebug(Info << "Enter getFaceCutPoints" << endl;)
    const faceList& faces = mesh_.faces();
    const pointField& points = mesh_.points();

    const labelList& pLabels = faces[fLabel];
    const label nPoints = pLabels.size();
    label pl1 = pLabels[0];
    scalar f1(f[0]);
//    if (mag(f1-f0) < 1e-12)
//    {
//        f1 = f0;
//    }
    forAll(pLabels,pi)
    {
        label pl2 = pLabels[(pi+1)%nPoints];
        scalar f2(f[(pi+1)%nPoints]);
        if (mag(f2-f0) < 1e-12)
        {
            f2 = f0;
        }
        bool edgeIsCut = min(f1,f2) < f0 && max(f1,f2) > f0;
        if ( edgeIsCut )
        {
            scalar s = (f0-f1)/(f2-f1);
            point pCut = points[pl1] + s*(points[pl2]-points[pl1]);
            cutPoints.append(pCut);
        }
        else if ( f1 == f0 )
        {
            cutPoints.append(points[pl1]);
        }
        pl1 = pl2;
        f1 = f2;
    }
}


void Foam::isoCutter::getFaceCutPoints
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
//    if (mag(f1-f0) < 1e-12)
//    {
//        f1 = f0;
//    }
    forAll(pts,pi)
    {
        label pi2 = (pi+1)%nPoints;
        scalar f2 = f[pi2];
        if (mag(f2-f0) < 1e-12)
        {
            f2 = f0;
        }
        bool edgeIsCut = min(f1,f2) < f0 && max(f1,f2) > f0;
        if ( edgeIsCut )
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

void Foam::isoCutter::fullySubmergedFaces
(
    const label cellI,
    const scalar f0,
    DynamicList<label>& fullSubFaces
)
{
//    const cellList& cells = mesh_.cells(); //It is important for the efficiency of the code to verify that this call does not generate the cell faces for all cells each time it is called but only the first time!
//    labelList cellFacesTest = cells[cellI];

//    const faceList& faces = mesh_.faces();
//    const labelList faceLabels = mesh_.cells()[cellI];
    const cell& fLabels = mesh_.cells()[cellI];

//  Info << "Into fullySubmergedFaces for cell " << cellI << " which has fLabels = " << fLabels << endl;
    forAll(fLabels,fi)
    {
        label fLabel = fLabels[fi];
        labelList pLabels = mesh_.faces()[fLabel];
//      Info << "fLabel = " << fLabel << ", pLabels = " << pLabels << endl;
        bool submerged(true);
        forAll(pLabels,pi)
        {
            label pLabel = pLabels[pi];
            if ( f_[pLabel] <= f0 ) //By convention all face points must be strictly below surface for the face to be fully submerged
            {
                submerged = false;
            }
        }
        if (submerged)
        {
            fullSubFaces.append(fLabel);
        }
    }
    fullSubFaces.shrink();
//  Info << "Out of fullySubmergedFaces with fullSubFaces = " << fullSubFaces << endl;
}


void Foam::isoCutter::writeFacesToPlyFile
(
    const List<pointField>& pfl,
    const word& fileName,
    const word& fileDir
)
{
    autoPtr<OFstream> plyFilePtr;
    plyFilePtr.reset(new OFstream(fileDir + "/" + fileName + ".ply"));
    plyFilePtr() << "ply" << endl;
    plyFilePtr() << "format ascii 1.0" << endl;
    plyFilePtr() << "comment " << fileName << endl;
    label nPoints(0);
    forAll(pfl,fi)
    {
        nPoints += pfl[fi].size();
    }

    plyFilePtr() << "element vertex " << nPoints << endl;
    plyFilePtr() << "property float32 x" << endl;
    plyFilePtr() << "property float32 y" << endl;
    plyFilePtr() << "property float32 z" << endl;
    plyFilePtr() << "element face " << pfl.size() << endl;
    plyFilePtr() << "property list uint8 int32 vertex_index" << endl;
    plyFilePtr() << "end_header" << endl;
    forAll(pfl,fi)
    {
        pointField pf = pfl[fi];
        forAll(pf,pi)
        {
            point p = pf[pi];
            plyFilePtr() << p[0] << " " << p[1] << " " << p[2] << endl;
        }
    }
    label np = 0;
    forAll(pfl,fi)
    {
        nPoints = pfl[fi].size();
        plyFilePtr() << nPoints;
        for (label pi = np; pi < np + nPoints; pi++ )
        {
            plyFilePtr() << " " << pi;
        }
        plyFilePtr() << "" << endl;
        np += nPoints;
    }
}


void Foam::isoCutter::subFaceFractions
(
    const scalar f0,
    surfaceScalarField& alphaf
)
{
    forAll(alphaf,fi)
    {
        subFaceFraction(fi,f0,alphaf[fi]);
    }
}


void Foam::isoCutter::subFaceFraction
(
    const label fi,
    const scalar f0,
    scalar& alphaf
)
{
    alphaf = 0;
    DynamicList<point> partSubFacePts((mesh_.faces()[fi]).size());
    bool fullySubmerged = getSubFace(fi,f0,partSubFacePts);
    if (fullySubmerged)
    {
        alphaf = 1;
    }
    else if (partSubFacePts.size() != 0)
    {
        vector fCtr, fArea;
        makeFaceCentreAndArea(partSubFacePts, fCtr, fArea);
        alphaf = mag(fArea)/mag(mesh_.Sf()[fi]);
    }
}

void Foam::isoCutter::isoFaceCentreAndArea
(
    const label ci,
    const scalar f0,
    vector& faceCentre,
    vector& faceArea
)
{
    DynamicList<label> cutFaces((mesh_.cells()[ci]).size());
    DynamicList<label> cutEdges(cutFaces.capacity());
    DynamicList<scalar> cutPoints(cutFaces.capacity());
    bool fullySubmerged(false);
    isoCutCell(ci,f0,cutFaces,cutEdges,cutPoints,fullySubmerged);
    DynamicList<point> isoPoints(cutFaces.size());
    getIsoFace(cutFaces,cutEdges,cutPoints,isoPoints);
    if (isoPoints.size()>0)
    {
        makeFaceCentreAndArea(isoPoints, faceCentre, faceArea);
    }
    else
    {
        Info << "Warning: isoPoints.size() <= 0 for cell " << ci << endl;
    }
}

void Foam::isoCutter::subCellFractions
(
    const scalar f0,
    volScalarField& alpha1
)
{
    for ( label ci = 0; ci < mesh_.nCells(); ci++ )
    {
        scalar ai;
        subCellFraction(ci, f0, ai);
        alpha1[ci] = ai;
    }
}

void Foam::isoCutter::vofCutCells
(
    const volScalarField& alpha1,
    const scalar tol,
    const label maxIter,
    volScalarField& f0
)
{
    vector subCellCtr;
    for ( label ci = 0; ci < mesh_.nCells(); ci++ )
    {
        scalar fi;
        scalar ai = alpha1[ci];
        if ( 0 < ai && ai < 1)
        {
            vofCutCell(ci, ai, tol, maxIter, fi, subCellCtr);
            f0[ci] = fi;
        }
        else
        {
            f0[ci] = ai;
        }
    }
}


void Foam::isoCutter::subCellFraction
(
    const label& ci,
    const scalar& f0,
    scalar& alpha1
)
{
    vector cellCtr;
    subCellFraction(ci, f0, alpha1, cellCtr);
}

void Foam::isoCutter::subCellFraction
(
    const label& ci,
    const scalar& f0,
    scalar& alpha1,
    vector& cellCtr
)
{
    alpha1 = 0;

	DynamicList<label> cutFaces((mesh_.cells()[ci]).size());
    DynamicList<label> cutEdges(cutFaces.capacity());
    DynamicList<scalar> cutPoints(cutFaces.capacity());

    bool fullySubmerged(false);
    isoCutCell(ci,f0,cutFaces,cutEdges,cutPoints,fullySubmerged);

    if ( cutFaces.size() > 0 )
    {
        //Get isoFace points
        DynamicList<point> isoPoints(cutFaces.size());
        getIsoFace(cutFaces,cutEdges,cutPoints,isoPoints);

        //Get points of cut cell face
        DynamicList< DynamicList<point> > partSubFacePoints(cutFaces.size());
        forAll(cutFaces,fi)
        {
			const label fLabel = cutFaces[fi];
            DynamicList<point> subFacePoints((mesh_.faces()[fLabel]).size());
            getSubFace(fLabel,f0,subFacePoints);
            partSubFacePoints.append(subFacePoints);
        }

        //Get points of fully submerged cell face
        DynamicList<label> fullSubFaces((mesh_.cells()[ci]).size());
        fullySubmergedFaces(ci,f0,fullSubFaces);

        DynamicList< DynamicList<point> > fullSubFacePoints(fullSubFaces.size());
        getFacePoints(fullSubFaces,fullSubFacePoints);

        //Gathering all cut cell face points in one list
        DynamicList< DynamicList<point> > cutCellFacePoints(partSubFacePoints);
        cutCellFacePoints.append(isoPoints);
        cutCellFacePoints.append(fullSubFacePoints);

        //Calculating cut cell face centres and areas
        vectorField faceCentres, faceAreas;
        makeFaceCentresAndAreas(cutCellFacePoints, faceCentres, faceAreas);

        //Calculating cut cell centre and volume
//      vector cellCtr;
        scalar cellVol;
        makeCellCentreAndVol(faceCentres, faceAreas, cellCtr, cellVol);
        const scalarField& V = mesh_.V();
        alpha1 = (cellVol/V[ci]);

    } else if ( fullySubmerged )
    {
        alpha1 = 1;
        cellCtr = mesh_.C()[ci];
    }
}

void Foam::isoCutter::makeFaceCentresAndAreas
(
    const List< DynamicList<point> >& pfList,
    vectorField& fCtrs,
    vectorField& fAreas
)
{
    DynamicList<vector> fCtrList(pfList.size()), fAreaList(pfList.size());
    forAll(pfList,pi)
    {
        vector fCtr, fArea;
        makeFaceCentreAndArea(pfList[pi], fCtr, fArea);
        fCtrList.append(fCtr);
        fAreaList.append(fArea);
    }
    fCtrs = fCtrList;
    fAreas = fAreaList;
}


void Foam::isoCutter::makeFaceCentreAndArea
(
    const DynamicList<point>& p,
    vector& fCtr,
    vector& fArea
)
{
    label nPoints = p.size();
    // If the face is a triangle, do a direct calculation for efficiency
    // and to avoid round-off error-related problems
    if ( nPoints == 3 )
    {
        fCtr = (1.0/3.0)*(p[0] + p[1] + p[2]);
        fArea = 0.5*((p[1] - p[0])^(p[2] - p[0]));
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
            fCtr = fCentre;
            fArea = vector::zero;
        }
        else
        {
            fCtr = (1.0/3.0)*sumAc/sumA;
            fArea = 0.5*sumN;
        }
    }
}


void Foam::isoCutter::getFaceCentresAndAreas
(
    const labelList& faceLabels,
    vectorField& fCtrs,
    vectorField& fAreas
)
{
    DynamicList<vector> fCtrList(faceLabels.size()), fAreaList(faceLabels.size());
    forAll(faceLabels,fi)
    {
        fCtrList.append(mesh_.Cf()[faceLabels[fi]]);
        fAreaList.append(mesh_.Sf()[faceLabels[fi]]);
    }
    fCtrs = fCtrList;
    fAreas = fAreaList;
}


void Foam::isoCutter::makeCellCentreAndVol
(
    const vectorField& fCtrs,
    const vectorField& fAreas,
    vector& cellCtr,
    scalar& cellVol
)
{
    // Clear the fields for accumulation
    cellCtr = vector::zero;
    cellVol = 0.0;

    // first estimate the approximate cell centre as the average of
    // face centres

    vector cEst(vector::zero);
    label nCellFaces(fCtrs.size());

    forAll(fCtrs, facei)
    {
        cEst += fCtrs[facei];
    }

    cEst /= nCellFaces;

    forAll(fCtrs, facei)
    {
        // Calculate 3*face-pyramid volume
        scalar pyr3Vol =
            max(mag(fAreas[facei] & (fCtrs[facei] - cEst)), VSMALL);

        // Calculate face-pyramid centre
        vector pc = (3.0/4.0)*fCtrs[facei] + (1.0/4.0)*cEst;

        // Accumulate volume-weighted face-pyramid centre
        cellCtr += pyr3Vol*pc;

        // Accumulate face-pyramid volume
        cellVol += pyr3Vol;
    }

    cellCtr /= cellVol;
    cellVol *= (1.0/3.0);
}


void Foam::isoCutter::printPoints
(
    const List<pointField>& pfl
)
{
    forAll(pfl,pfli)
    {
        Info << "p{" << pfli+1 << "} = [";
        printPoints(pfl[pfli]);
        Info << "];" << endl;
    }
}


void Foam::isoCutter::printPoints
(
    const pointField& pf
)
{
    forAll(pf,pfi)
    {
        point p = pf[pfi];
        Info << "[";
        forAll(p,pi)
        {
            Info << " " << p[pi];
        }
        Info << "];";
    }
}


void Foam::isoCutter::getFacePoints
(
    const labelList& fLabels,
    DynamicList< DynamicList<point> >& facePointLists
)
{
    forAll(fLabels,fi)
    {
		const label fLabel = fLabels[fi];
        DynamicList<point> fp((mesh_.faces()[fLabel]).size());
        getFacePoints(fLabel,fp);
        facePointLists.append(fp);
    }
}


void Foam::isoCutter::getFacePoints
(
    const label faceLabel,
    DynamicList<point>& fp
)
{
    const pointField& points = mesh_.points();
    const labelList& fpLabels = mesh_.faces()[faceLabel];
    forAll(fpLabels,pi)
    {
        fp.append(points[fpLabels[pi]]);
    }
}

// ************************************************************************* //
