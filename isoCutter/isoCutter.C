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
#include "volPointInterpolation.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isoCutter::isoCutter
(
    const fvMesh& mesh,
    const scalarField& f,
    const scalar& f0
)
:
    mesh_(mesh),
    f_(f),
    f0_(f0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::isoCutter::vofCutCell
(
	const label& ci,
	const scalar& alpha1,
	const scalar& tol,
	const label& maxIter,
	scalar& f0,
	vector& subCellCtr
)
{
	Info << "Enter vofCutCell for cell " << ci << " which has alpha1 = " << alpha1 << endl;
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
//	Info << "fMin = " << fMin << " fMax = " << fMax << endl;
	
	//Initial guess of isovalue
	scalar aMin(0), aMax(1), alpha0;
	f0 = (alpha1 - aMin)/(aMax-aMin)*(fMin-fMax) + fMax;
//	f0 = 0.5*(fMin + fMax);
	subCellFraction(ci, f0, alpha0, subCellCtr);
//	Info << "f0 = " << f0 << " gives alpha = " << alpha0 << endl;
	
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
//		f0 = (alpha1 - aMin)/(aMax-aMin)*(fMin-fMax) + fMax; //This does an extremely poor job in narrowing in the interval - especially for almost filled cells
		subCellFraction(ci, f0, alpha0, subCellCtr);
//		Info << nIter << ": f0 = " << f0 << " gives alpha = " << alpha0 << endl;
		nIter++;
	}
	Info << nIter-1 << ": f0 = " << f0 << " gives alpha = " << alpha0 << endl;
}

void Foam::isoCutter::isoCutCell
(
    const label& cellI,
	const scalar& f0,
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
    const label& cellI,
	const scalar& f0,
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
    const label& cellI,
	const scalar& f0,
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
    const label& cellI,
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
    pointField& isoPoints
)
{
    const faceList& faces = mesh_.faces();
    const pointField& points = mesh_.points();
    bool findDuplicates = false;

    pointField p;
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
        label nPoints = p.size();
        forAll(p,pi)
        {
            if ( mag(p[pi]-p[(pi+1)%nPoints]) > VSMALL )
            {
                isoPoints.append(p[pi]);
            }
        }
        if ( isoPoints.size() < 1 && p.size() > 0 )
        {
            isoPoints.append(p[0]);
        }
    }
    else
    {
        isoPoints = p;
    }
}

bool Foam::isoCutter::getSubFace
(
    const label& faceLabel,
	const scalar& f0,
    pointField& partSubFacePts
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

void Foam::isoCutter::fullySubmergedFaces
(
    const label& cellI,
	const scalar& f0,
    DynamicList<label>& fullSubFaces
)
{
    const faceList& faces = mesh_.faces();
    const labelList faceLabels = mesh_.cells()[cellI];

    forAll(faceLabels,fi)
    {
        labelList face = faces[faceLabels[fi]];
        bool submerged(true);
        forAll(face,pi)
        {
            if ( f_[face[pi]] <= f0 ) //By convention all face points must be strictly below surface for the face to be fully submerged
            {
                submerged = false;
            }
        }
        if (submerged)
        {
            fullSubFaces.append(faceLabels[fi]);
        }
    }
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
	const scalar& f0,
    surfaceScalarField& alphaf
)
{
    alphaf = 0;
    const vectorField& Sf = mesh_.Sf();

    forAll(Sf,fi)
    {
        pointField partSubFacePts;
        bool fullySubmerged = getSubFace(fi,f0,partSubFacePts);
        if (fullySubmerged)
        {
            alphaf[fi] = 1;
        }
        else if (partSubFacePts.size() != 0)
        {
            vector fCtr, fArea;
            makeFaceCentreAndArea(partSubFacePts, fCtr, fArea);
            alphaf[fi] = mag(fArea)/mag(Sf[fi]);
        }
    }
}


void Foam::isoCutter::updateAlpha
(
	const volScalarField& alpha,
	const surfaceScalarField& phi,
	const volVectorField& U,
	const scalar& dt,
	scalarField& dV
)
{
	Info << "Enter updateAlpha" << endl;
	volPointInterpolation vpi(mesh_);
	const scalarField alphap = vpi.interpolate(alpha);
    const cellList& cells = mesh_.cells();
	dV = 0.0; //estimated total water volume transported across mesh faces during time interval dt

	forAll(cells,ci)
	{	
		//Make list of all cell faces out of which fluid is flowing
		DynamicList<label> outFluxingFaces;
		getOutFluxFaces(phi, ci, outFluxingFaces);

		//Determine whether cell is fully above or fully below surface (if neither, it must be cut)
		scalar aMin, aMax;
		subSetExtrema(alphap,mesh_.cellPoints()[ci],aMin,aMax);

		//Set dV on all outfluxing faces in accordance to the status of the cell (if allPointsAbove dV should remain 0 on all outfluxing faces so do nothing)
//		if ( allPointsBelow || (alpha[ci] <= 0 || alpha[ci] >= 1) )
		if ( aMin > 0.5 )
		{
			forAll(outFluxingFaces,fi)
			{
				label lfi = outFluxingFaces[fi];
				dV[lfi] = phi[lfi]*dt;
			}
		}
//		else if ((!allPointsAbove) && ( 0<alpha[ci] && alpha[ci] < 1 )) //Then cell must be cut
		else if ( aMax > 0.5 ) //Then cell must be cut
		{
			Info << "------------ Cell " << ci << " ------------" << endl;
			Info << "outFluxingFaces: " << outFluxingFaces << endl;
			//Calculate isovalue f0 reproducing VOF value to given accuracy
			scalar f0, tol(1e-3);
			label maxIter(100);
			vector subCellCtr;
			vofCutCell(ci, alpha[ci], tol, maxIter, f0, subCellCtr);
			Info << "Isovalue giving proper VOF value: f0 " << f0 << ", subCellCtr: " << subCellCtr << endl;

			//Calculate isoFace0 centre xs0, normal ns0, and velocity U0 = U(xs0)
			vector x0(vector::zero), n0(vector::zero);
			isoFaceCentreAndArea(ci,f0,x0,n0);
			Info << "isoFace0 centre and area are x0 = " << x0 << " and n0 = " << n0 << endl; 
			//Make n0 point out of subCell
			if (((x0 - subCellCtr) & n0) < 0)
			{
				n0 *= (-1.0);
				Info << "Changing direction of n0 " << n0 << endl;
			}
			
			//Interpolate velocity to isoFace centre
			interpolationCellPoint<vector> UInterp(U);
			vector U0 = UInterp.interpolate(x0,ci);
			scalar Un0 = (U0 & n0);
			if (mag(n0) > 0)
			{
				Un0 /= mag(n0);
			}
			else
			{
				Info << "WARNING: n0 has zero length. Probably because isoFace is a point. Skipping division of Un0 by mag(n0)." << endl;
			}
			Info << "Velocity interpolated to isoFace0 centre, U0 = " << U0 << " and projected on n0, Un0 = " << Un0 << endl;
			
			vector x00(x0), n00(n0); //Necessary when more than one outfluxing face
			//Estimating time integrated water flux through each outfluxing face
			forAll(outFluxingFaces,fi)
			{
				label fLabel = outFluxingFaces[fi];
				Info << "FACE " << fLabel << " ------------" << endl;
	
//------------------------Calculate initial alphaf on the face------------------------
				scalar t0 = 0.0; //Initial time
				scalar a0 = 0.0; //Initial alphaf
				x0 = x00; //Necessary when more than one outfluxing face
				n0 = n00; //Necessary when more than one outfluxing face
				const vectorField& Sf = mesh_.Sf();
				pointField partSubFacePts;
				bool fullySubmerged = getSubFace(fLabel,f0,partSubFacePts);
				if (fullySubmerged)
				{
					a0 = 1;
				}
				else if (partSubFacePts.size() != 0)
				{
					vector fCtr, fArea;
					makeFaceCentreAndArea(partSubFacePts, fCtr, fArea);
					a0 = mag(fArea)/mag(Sf[fLabel]);
				}
				Info << "a0 on face: " << a0 << endl;
				
				//Generating list of face vertices approached by isoFace
				scalarField av; //Unique version of aVert
				vector xFirst(vector::zero), xLast(vector::zero);
				getVertexAlphas(alphap,fLabel,f0,Un0,av,xFirst,xLast);
//----------------------------------------------------------------------
				
				if (av.size() > 0) //((face is either completely submerged OR unsubmerged) AND Un is such that surface moves away from face) OR (face is cut and Un = 0).
				{
					label na = 0;
					scalar a1(0.0); //VOF value on face corresponding to isoFace1
					scalar t1(0.0); //Next time a vertex is met
					while (t0 < dt && na < av.size())
					{
						Info << "Calculating flux integral contribution for sub time step " << na << " starting at t0 " << t0 << " (dt = " << dt << ")" << endl;
						//Calculating new isoFace position and orientation
						scalar f1(av[na]);
						Info << "Calculating isoFace1 for f1 = " << f1 << endl;
						vector x1(vector::zero), n1(vector::zero); //IsoFace1 centre and area vectors
						//scalar a1(0.0); //VOF value on face corresponding to isoFace1
						
						if (na == 0 && (a0 == 0 || a0 == 1)) //If face is initially completely full or empty so will it be at its first encounter with surface
						{
							a1 = a0;
							x1 = xFirst;
							n1 = n0;
							Info << "Since na = " << na << " and a0 = " << a0 << " we set a1 = a0 and x1 = xFirst = " << xFirst << endl;
						}
						else if (na == av.size()-1) //At the last point the face will either be completely filled or emptied - so isoFace = the vertex point
						{
							a1 = pos(f0-f1); //If f0 > f1, face is filling up with water - else it is being emptied
							x1 = xLast;
							n1 = n0; //since isoFace1 is a point its normal is undefined and we set it to previous value
							Info << "Since na = " << na << " = na == av.size()-1 we set a1 = pos(f0-f1) = " << pos(f0-f1) << " and x1 = xLast = " << xLast << endl;
						} 
						else
						{
							isoFaceCentreAndArea(ci,f1,x1,n1);
	
							//Finding alphaf on face for new time
							fullySubmerged = getSubFace(fi,f0,partSubFacePts);
							if (fullySubmerged)
							{
								a1 = 1;
							}
							else if (partSubFacePts.size() != 0)
							{
								vector fCtr, fArea;
								makeFaceCentreAndArea(partSubFacePts, fCtr, fArea);
								a1 = mag(fArea)/mag(Sf[fLabel]);
							}
						}
						Info << "isoFace1 centre and area are x1 = " << x1 << " and n1 = " << n1 << " has a1 = " << a1 << endl; 
/*						if ((n1 & n0) < 0)
						{
							n1 *= (-1.0);
							Info << "Changing direction of n1: " << n1 << endl;
						}
						
						//Interpolate velocity to isoFace centre - currently not used
						vector U1 = UInterp.interpolate(x1,ci);
						scalar Un1 = (U1 & n1)/mag(n1);
						Info << "Velocity interpolated to isoFace1 centre, U1 = " << U1 << " and projected on n1, Un1 = " << Un1 << endl;
*/						
						//Estimating time where the surface will be at x1
						scalar dtn = ((x1-x0) & (n0/mag(n0)))/Un0; //What about when Un0 = 0?
						Info << "x0: " << x0 << " x1: " << x1 << " n0: " << n0 << " U0: " << U0 << endl;
						t1 = min(t0 + dtn,dt);
						Info << "Estimated dtn " << dtn << " and t1 " << t1 << endl;

						a1 = a0 + (t1-t0)/dtn*(a1-a0); //Linear interpolation of alphaf to time t1 from values at time t0 and t0+dtn - only changes a1 if t0+dtn > dt.
						Info << "Resulting a1: " << a1 << endl;
						dV[fLabel] += phi[fLabel]*0.5*(a0 + a1)*(t1-t0); //Trapezoidal estimate
						Info << "dV[fLabel] = " << dV[fLabel] << " after adding " << phi[fLabel]*0.5*(a0 + a1)*(t1-t0) << " m3" << endl;
						
						t0 = t1;
						x0 = x1;
//						n0 = n1;
						a0 = a1;
//						U0 = U1;
//						Un0 = Un1;
						na++;
					}
					if (t1 < dt)
					{
						dV[fLabel] += phi[fLabel]*a1*(dt-t1);
					}
				}
				else
				{
					Info << "No elements in aVert. Face must be uncut and surface moving away from it." << endl;
					dV[fLabel] += phi[fLabel]*a0*dt;
				}

			}
		}
	}
}

void Foam::isoCutter::getOutFluxFaces
(
	const surfaceScalarField& phi,
	const label& ci,
	DynamicList<label>& outFluxingFaces
)
{
    const cellList& cells = mesh_.cells();
	forAll(cells[ci],fi)
	{
		label fLabel = cells[ci][fi];
		if (mesh_.owner()[fLabel] == ci && phi[fLabel] > 0)
		{
			outFluxingFaces.append(fLabel);
		}
		else if (ci < mesh_.nInternalFaces())
		{
			if (mesh_.neighbour()[fLabel] == ci && phi[fLabel] < 0)
			{
				outFluxingFaces.append(fLabel);					
			}
		}
	}
	outFluxingFaces.shrink();
}

void Foam::isoCutter::getVertexAlphas
(
	const scalarField& alphap,
	const label& fLabel,
	const scalar& f0,
	const scalar& Un0,
	scalarField& av,
	vector& xFirst,
	vector& xLast
)
{
	labelList pLabels = mesh_.faces()[fLabel];
//	Info << "pLabels: " << pLabels << endl;

	vectorField xVert; //Approached face vertex points
	scalarField aVert; //Corresponding alpha
	forAll(pLabels,pi)
	{
		scalar aVerti = alphap[pLabels[pi]];
//	Info << "aVerti for pi = " << pi << " equals " << aVerti << endl;
		if ( neg( Un0*(aVerti - f0) ) ) //Neg if Un0 and (aVerti-f0) have opposite signs i.e. if ( Un0 > 0 and (aVerti<f0) ) or ( Un0 < 0 and (aVerti>f0) )
		{
			xVert.append(mesh_.points()[pLabels[pi]]);
			aVert.append(aVerti);
		}
	}
	Info << "aVert = " << aVert << ", xVert = " << xVert << endl;
	
	if (aVert.size() > 0) //((face is either completely submerged OR unsubmerged) AND Un is such that surface moves away from face) OR (face is cut and Un = 0).
	{
		scalar aMin(GREAT), aMax(-GREAT);
		label iaMin(-1), iaMax(-1);
		forAll(aVert,ai)
		{
			if (aVert[ai] < aMin)
			{
				aMin = aVert[ai];
				iaMin = ai;
			}
			if (aVert[ai] > aMax)
			{
				aMax = aVert[ai];
				iaMax = ai;
			}
		}
		xFirst = xVert[iaMax];
		xLast = xVert[iaMin];
		if (Un0 < 0)
		{
			Swap(xFirst,xLast);
		}
		sort(aVert);
		
//					Info << "Sorted aVert: " << aVert << endl;
		av.append(aVert[0]);
		for (label n = 1; n < aVert.size(); n++)
		{
			Info << "mag(aVert[n] - aVert[n-1]): " << mag(aVert[n] - aVert[n-1]) << endl;
			if (mag(aVert[n] - aVert[n-1]) > 10*SMALL) //SMALL is 1e-15 and I have experienced the difference being 1.2 e-15.
			{
				av.append(aVert[n]);
			}
		}
//					Info << "Unique values of aVert in av: " << av << endl;
		if (Un0 > 0) //Water moves "forward", so nearest vertices are those with highest alphaf value
		{
			reverse(av);
		}
	}
	Info << "Reversed av: " << av << endl;	
}


void Foam::isoCutter::subSetExtrema
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

void Foam::isoCutter::isoFaceCentreAndArea
(
	const label& ci,
	const scalar& f0,
	vector& faceCentre,
	vector& faceArea
)
{
	DynamicList<label> cutFaces;
	DynamicList<label> cutEdges;
	DynamicList<scalar> cutPoints;
	bool fullySubmerged(false);
	isoCutCell(ci,f0,cutFaces,cutEdges,cutPoints,fullySubmerged);
	pointField isoPoints;
	getIsoFace(cutFaces,cutEdges,cutPoints,isoPoints);
	if (isoPoints.size()>0)
	{
		makeFaceCentreAndArea(isoPoints, faceCentre, faceArea);
	}
	else
	{
		Info << "Warning: isoPoints.size() <= 0" << endl;
	}	
}

void Foam::isoCutter::subCellFractions
(
	const scalar& f0,
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
	const scalar& tol,
	const label& maxIter,
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

	DynamicList<label> cutFaces;
	DynamicList<label> cutEdges;
	DynamicList<scalar> cutPoints;
	DynamicList<labelList> submergedPoints;
	bool fullySubmerged(false);
	isoCutCell(ci,f0,cutFaces,cutEdges,cutPoints,fullySubmerged);

	if ( cutFaces.size() > 0 )
	{
		//Get isoFace points
		pointField isoPoints;
		getIsoFace(cutFaces,cutEdges,cutPoints,isoPoints);

		//Get points of cut cell face
		DynamicList<pointField> partSubFacePoints;
		forAll(cutFaces,fi)
		{
			pointField subFacePoints;
			getSubFace(cutFaces[fi],f0,subFacePoints);
			partSubFacePoints.append(subFacePoints);
		}

		//Get points of fully submerged cell face			
		DynamicList<label> fullSubFaces;
		fullySubmergedFaces(ci,f0,fullSubFaces);
		DynamicList<pointField> fullSubFacePoints;
		getFacePoints(fullSubFaces,fullSubFacePoints);
		
		//Gathering all cut cell face points in one list
		DynamicList<pointField> cutCellFacePoints(partSubFacePoints);
		cutCellFacePoints.append(isoPoints);
		cutCellFacePoints.append(fullSubFacePoints);

		//Calculating cut cell face centres and areas
		vectorField faceCentres, faceAreas;
		makeFaceCentresAndAreas(cutCellFacePoints, faceCentres, faceAreas);

		//Calculating cut cell centre and volume
//		vector cellCtr;
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


void Foam::isoCutter::write()
{
    const vectorField& C = mesh_.C();
    const scalarField& V = mesh_.V();
    const vectorField& Cf = mesh_.Cf();
    const vectorField& Sf = mesh_.Sf();

    DynamicList<pointField> isoFaces;
    DynamicList<pointField> pSubFaces;
    DynamicList<pointField> subFaces;

    for ( label ci = 0; ci < mesh_.nCells(); ci++ )
    {
        DynamicList<label> cutFaces;
        DynamicList<label> cutEdges;
        DynamicList<scalar> cutPoints;
        bool fullySubmerged(false);
        isoCutCell(ci,f0_,cutFaces,cutEdges,cutPoints,fullySubmerged); //NOTE THAT THIS USES F0_ SO NOT NECESSARILY PRINTING THE RIGHT THING!!!

        if ( cutFaces.size() > 0 )
        {
//            Info << "Cell: " << ci << endl;

			//Get isoFace points
            pointField isoPoints;
            getIsoFace(cutFaces,cutEdges,cutPoints,isoPoints);

			//Get points of cut cell face
            DynamicList<pointField> partSubFacePoints;
            forAll(cutFaces,fi)
            {
                pointField subFacePoints;
                getSubFace(cutFaces[fi],f0_, subFacePoints);
                partSubFacePoints.append(subFacePoints);
            }

			//Get points of fully submerged cell face			
            DynamicList<label> fullSubFaces;
            fullySubmergedFaces(ci,f0_,fullSubFaces);
			DynamicList<pointField> fullSubFacePoints;
			getFacePoints(fullSubFaces,fullSubFacePoints);
			
			//Gathering all cut cell face points in one list
			DynamicList<pointField> cutCellFacePoints(partSubFacePoints);
			cutCellFacePoints.append(isoPoints);
			cutCellFacePoints.append(fullSubFacePoints);

			//Calculating cut cell face centres and areas
            vectorField faceCentres, faceAreas;
            makeFaceCentresAndAreas(cutCellFacePoints, faceCentres, faceAreas);

			//Calculating cut cell centre and volume
            vector cellCtr;
            scalar cellVol;
            makeCellCentreAndVol(faceCentres, faceAreas, cellCtr, cellVol);

			//Printing cut faces
			bool printCellAndFaceInfoToLog = false;
			if ( printCellAndFaceInfoToLog )
			{
				Info << "C = " << C[ci] << ", C_cut " << cellCtr << ", V = " << V[ci] << ", V_cut = " << cellVol << ", V_cut/V = "  << (cellVol/V[ci]) << endl;
				forAll(cutFaces,fi)
				{
					Info << "Cut face: " << cutFaces[fi] 
						 << ", Sf = " << Sf[cutFaces[fi]] << ", Sf_cut = " << faceAreas[fi] 
						 << ", Cf = " << Cf[cutFaces[fi]] << ", Cf_cut = " << faceCentres[fi] 
						 << " |Sf_cut|/|Sf| = " << mag(faceAreas[fi])/mag(Sf[cutFaces[fi]]) 
						 << endl;
				}
				Info << "Cut cell faces points (first " << cutFaces.size() 
					 << " are the cut faces, the next is the isoFace, and any remaining are fully submerged faces):" 
					 << endl;
				printPoints(cutCellFacePoints);
				Info << "---------------------------" << endl;
			}
			
			//Writing submerged cell faces to ply files
			isoFaces.append(isoPoints);
			pSubFaces.append(partSubFacePoints);
			subFaces.append(fullSubFacePoints);

        }
    }
    writeFacesToPlyFile(isoFaces, "isoFaces", mesh_.time().timeName());
    writeFacesToPlyFile(pSubFaces, "cutFaces", mesh_.time().timeName());
    writeFacesToPlyFile(subFaces, "subFaces", mesh_.time().timeName());
}


void Foam::isoCutter::makeFaceCentresAndAreas
(
    const List<pointField>& pfList,
    vectorField& fCtrs,
    vectorField& fAreas
)
{
    forAll(pfList,pi)
    {
        vector fCtr, fArea;
        makeFaceCentreAndArea(pfList[pi], fCtr, fArea);
        fCtrs.append(fCtr);
        fAreas.append(fArea);
    }
}


void Foam::isoCutter::makeFaceCentreAndArea
(
    const pointField& p,
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
    forAll(faceLabels,fi)
    {
        fCtrs.append(mesh_.Cf()[faceLabels[fi]]);
        fAreas.append(mesh_.Sf()[faceLabels[fi]]);
    }
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
	DynamicList<pointField>& facePointLists
)
{
	forAll(fLabels,fi)
	{
		pointField fp;
		getFacePoints(fLabels[fi],fp);
		facePointLists.append(fp);
	}
}


void Foam::isoCutter::getFacePoints
(
	const label faceLabel,
	pointField& fp
)
{
	const pointField& points = mesh_.points();

	labelList fpLabels = mesh_.faces()[faceLabel];
	forAll(fpLabels,pi)
	{
		fp.append(points[fpLabels[pi]]);
	}	
}

// ************************************************************************* //
