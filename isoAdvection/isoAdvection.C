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

	//Upwinding on all internal faces receiving fluid from a non-surface cell
	forAll(dVf,fi)
	{
		label donorCell(-1);
		if ( phi_[fi] >= 0 )
		{
			donorCell = mesh_.owner()[fi];
		}
		else if ( fi < mesh_.nInternalFaces() )
		{
			donorCell = mesh_.neighbour()[fi];
		}
		if (donorCell != -1)
		{
			dVf[fi] = phi_[fi]*alpha1_[donorCell]*dt;
		}
		else
		{
			Info << "Warning: No donor cell for face " << fi << endl;
		}
	}
	
	//Isocutting estimate of water flux from surface cellss
    forAll(surfaceCells,cellI)
    {
        const label ci = surfaceCells[cellI];
        Info << "\n------------ Cell " << ci << " with alpha1_ = " << alpha1_[ci] << " ------------" << endl;

        //Make list of all cell faces out of which fluid is flowing
        DynamicList<label> outFluxingFaces;
		getOutFluxFaces(ci, outFluxingFaces);
        Info << "outFluxingFaces: " << outFluxingFaces << "\n" << endl;

        //Calculate isoFace0 centre xs0, normal ns0, and velocity U0 = U(xs0)
        scalar f0(0.0), Un0(0.0);
        vector x0(vector::zero), n0(vector::zero);
        calcIsoFace(ci,x0,n0,f0,Un0); //This one really also should give us a0 on all faces since it is calculated anyway. Do this with a cutCell structure
		Info << "calcIsoFace gives initial surface: \nx0 = " << x0 << ", \nn0 = " << n0 << ", \nf0 = " << f0 << ", \nUn0 = " << Un0 << endl; 

        //Estimating time integrated water flux through each outfluxing face
        forAll(outFluxingFaces,fi)
        {
            const label fLabel = outFluxingFaces[fi];
			Info << "\nFace " << fLabel << " with outward normal nf = " << sign( (mesh_.Cf()[fLabel]-mesh_.C()[ci]) & mesh_.Sf()[fLabel] )*mesh_.Sf()[fLabel]/mag(mesh_.Sf()[fLabel]) << endl;
            dVf[fLabel] = timeIntegratedFlux(fLabel,x0,n0,Un0,f0,dt);
        }
    }
}

void Foam::isoAdvection::findSurfaceCells
(
    DynamicList<label>& surfaceCells
)
{
    isSurfaceCell_ = false;
    forAll(alpha1_,ci)
    {
        scalar aMin(GREAT), aMax(-GREAT);
        subSetExtrema(ap_,mesh_.cellPoints()[ci],aMin,aMax);
        if ( (aMin < 0.5 && aMax > 0.5) )
//        if ( (aMin < 0.5 && aMax > 0.5) && ( 1e-6 < alpha1_[ci] && alpha1_[ci] < 1-1e-6 ) )
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
	Info << "Enter calcIsoFace" << endl;
    //Construction of isosurface calculator to get access to its functionality
    isoCutter cutter(mesh_,ap_);
    scalar tol(1e-6);
    label maxIter(100);
    vector subCellCtr;
    cutter.vofCutCell(ci, alpha1_[ci], tol, maxIter, f0, subCellCtr);

		Info << "Cell is neither almost full or almost empty" << endl;
        cutter.isoFaceCentreAndArea(ci,f0,x0,n0); //Stupid to recalculate this here - should be provided by vofCutCell above
		
		if (mag(n0)<SMALL) //Cell almost full or empty so isoFace is undefined. Calculating normal by going a little into the cell
		{
			scalar aMin(GREAT), aMax(-GREAT);
			subSetExtrema(ap_,mesh_.cellPoints()[ci],aMin,aMax);
			if (mag(f0-aMin) < mag(f0-aMax)) //f0 almost equal to aMin, i.e. cell almost full
			{
				Info << "Cell is almost full" << endl;
				scalar f0Inside =  aMin + tol*(aMax-aMin);
				cutter.isoFaceCentreAndArea(ci,f0Inside,x0,n0);
				if (((x0 - mesh_.C()[ci]) & n0) < 0.0) //n0 should point from cell centre towards isoFace centre for an almost full cell
				{
					n0 *= (-1.0);
					Info << "Changing direction of n0 " << n0 << endl;
				}
			}
			else //f0 almost equal to aMax, i.e. cell almost empty
			{
				Info << "Cell is almost empty" << endl;
				scalar f0Inside =  aMax + tol*(aMin-aMax);
				cutter.isoFaceCentreAndArea(ci,f0Inside,x0,n0);
				if (((mesh_.C()[ci] - x0) & n0) < 0.0) //n0 should point from isoFace centre towards cell centre for an almost empty cell
				{
					n0 *= (-1.0);
					Info << "Changing direction of n0 " << n0 << endl;
				}
			}

		}
        //Make n0 point out of subCell
        else if (((x0 - subCellCtr) & n0) < 0) //n0 should point out of water surface, i.e. in the direction from subCell centre to isoFace centre
        {
            n0 *= (-1.0);
            Info << "Changing direction of n0 " << n0 << endl;
        }

		Info << "Normalising n0: " << n0 << endl;
	n0 /= mag(n0);

    //Interpolate velocity to isoFace centre
    interpolationCellPoint<vector> UInterp(U_);
    vector U0 = UInterp.interpolate(x0,ci);
    Un0 = (U0 & n0);
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
	Info << "sortedTimes = " << sortedTimes << endl;
	
    scalar dVf(0.0); //Volume flowing through face in time interval [0,dt] to be calculated below

    //Dealing with case where face is not cut by surface during time interval [0,dt] because face was already passed by surface
    if ( sortedTimes[nPoints-1] <= 0.0 ) //All cuttings in the past
    {
		Info << "All cuttings in the past" << endl;
        dVf = phi_[fLabel]*dt*pos(Un0); //If all face cuttings were in the past and cell is filling up (Un0>0) then face must be full during whole time interval
        return dVf;
    }

    //Dealing with case where face is not cut by surface during time interval [0,dt] because dt is too small for surface to reach closest face point
    if ( sortedTimes[0] >= dt ) //All cuttings in the future
    {
		Info << "All cuttings in the future" << endl;
        dVf = phi_[fLabel]*dt*(1-pos(Un0)); //If all cuttings are in the future but non of them within [0,dt] then if cell is filling up (Un0 > 0) face must be empty during whole time interval
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
	Info << "Cutting sortedTimes at 0 and dt: t = " << t << endl;	
	
    bool faceUncutInFirstInterval(sortedTimes[0] > 0.0);
    bool faceUncutInLastInterval(sortedTimes[nPoints-1] < dt);

    //Dealing with cases where face is cut at least once during time interval [0,dt]
    label nt = 0; //Time interval counter
    scalar initialArea(0.0); //Submerged area at time
    if ( faceUncutInFirstInterval ) //Special treatment for first time interval if face is uncut during this
    {
        dVf = phi_[fLabel]*(t[nt+1]-t[nt])*(1.0-pos(Un0)); //If Un0 > 0 cell is filling up - hence if face is cut at a later time but not initially it must be initially empty
        initialArea = mag(mesh_.Sf()[fLabel])*(1.0-pos(Un0));
		Info << "faceUncutInFirstInterval, so special treatment for first time interval: [" << t[nt] << ", " << t[nt+1] << "] giving dVf = " << dVf << endl;
        nt++;
    }
	else //calculate initialArea if face is initially cut
	{
		isoCutter cutter(mesh_,ap_);
		initialArea = cutter.getSubFaceFraction(fLabel, -sign(Un0)*pTimes, 0.0)*mag(mesh_.Sf()[fLabel]); //does not use ap_
        scalar A(0.0), B(0.0);
        quadAreaCoeffs(fLabel,pTimes,t[nt],t[nt+1],A,B);
	}

	Info << "InitialArea for next time step corresponds to face phase fraction a0 = " << initialArea/mag(mesh_.Sf()[fLabel]) << " where |Sf| = " << mag(mesh_.Sf()[fLabel]) << " was used." << endl;
    while ( nt < t.size()-(1+faceUncutInLastInterval) ) //
    {
        scalar A(0.0), B(0.0);
        quadAreaCoeffs(fLabel,pTimes,t[nt],t[nt+1],A,B);
        scalar integratedQuadArea = sign(Un0)*(A/3.0 + 0.5*B); //Integration of area(t) = A*t^2+B*t from t = 0 to 1.
        scalar Unf = phi_[fLabel]/mag(mesh_.Sf()[fLabel]); //normal velocity on face
        dVf += Unf*(t[nt+1]-t[nt])*(initialArea + integratedQuadArea);
        initialArea += sign(Un0)*(A + B); //Adding quad area to submerged area
		Info << "Integrating area for " << nt+1 << "'th time interval: [" << t[nt] << ", " << t[nt+1] << "] giving dVf = " << dVf << " and a0 = " << initialArea/mag(mesh_.Sf()[fLabel]) << endl;
		Info << "face owner = " << mesh_.owner()[fLabel] << endl;
//      scalar newdVf = phi_[fLabel]/mag(mesh_.Sf()[fLabel])*(t1-t0)*(a0*mag(mesh_.Sf()[fLabel]) + sign(Un0)*integratedArea(fLabel,f0,f1));
        nt++;
    }

    if ( faceUncutInLastInterval ) //Special treatment for last time interval if face is uncut during this
    {
		Info << "faceUncutInLastInterval, so special treatment for last (" << nt+1 << "'th) time interval: [" << t[nt] << ", " << t[nt+1] << "]" << endl;
        dVf += phi_[fLabel]*(t[nt+1]-t[nt])*pos(Un0); //If face is cut at some intermediate time but not at last time, then if Un0 > 0 (cell filling up) face must be filled at last time interval.
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
	const labelList fLabels = cells[ci];
    forAll(fLabels,fi)
    {
        const label fLabel = fLabels[fi];
		if (fLabel < mesh_.nInternalFaces())
		{
			const label owner = mesh_.owner()[fLabel];
			if (owner == ci)
			{
				if (phi_[fLabel] > 0.0)
				{
					outFluxingFaces.append(fLabel);
				}
			}
			else if ( phi_[fLabel] < 0.0 ) //ci must be neighbour of fLabel
			{
				outFluxingFaces.append(fLabel);
			}
		}
    }
    outFluxingFaces.shrink();
}



Foam::label Foam::isoAdvection::otherCell
(
	const label& fLabel,
	const label& cLabel
)
{
	label otherCellLabel(-1);
	if (fLabel < mesh_.nInternalFaces())
	{
		if (cLabel == mesh_.owner()[fLabel])
		{
			otherCellLabel = mesh_.neighbour()[fLabel];
		}
		else
		{
			otherCellLabel = mesh_.owner()[fLabel];
		}
	}
	return otherCellLabel;
}

/*
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
*/

void Foam::isoAdvection::quadAreaCoeffs
(
    const label& fLabel,
    const scalarField& f,
    const scalar& f0,
    const scalar& f1,
    scalar& alpha,
    scalar& beta
)
{
    isoCutter cutter(mesh_,ap_);
    pointField pf0, pf1;
    cutter.getFaceCutPoints(fLabel,f,f0,pf0);
    cutter.getFaceCutPoints(fLabel,f,f1,pf1);
//  cutter.getFaceCutPoints(fLabel,f0,pf0);
//  cutter.getFaceCutPoints(fLabel,f1,pf1);

    label np0(pf0.size()), np1(pf1.size());
//    Info << "Face " << fLabel << " was cut " << np0 << " times by f0 = " << f0 << " and " << np1 << " times by " << " f1 = " << f1 << endl;

//  scalar area(0.0);
    alpha = 0.0;
    beta = 0.0;

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


		//Defining quadrilateral vertices
		vector A(pf0[0]), C(pf1[0]), B(vector::zero), D(vector::zero);

		//Triangle cases
        if (np0 == 2)
		{
			B = pf0[1];
		}
		else
		{
			B = A;
		}

        if (np1 == 2)
		{
			D = pf1[1];
		}
		else
		{
			D = C;
		}
/*
        if (np0 == 1)
        {
			Info << "Warning: np0 == 1" << endl;
			Info << "xhat = " << xhat << ", yhat = " << yhat << ", zhat = " << zhat << ". x.x = " << (xhat & xhat) << ", y.y = " << (yhat & yhat) <<", z.z = " << (zhat & zhat) << ", x.y = " << (xhat & yhat) << ", x.z = " << (xhat & zhat) << ", y.z = " << (yhat & zhat) << endl;
            pf0.append(pf0[0]);
        }
        else if (np1 == 1)
        {
			Info << "Warning: np1 == 1" << endl;
			Info << "xhat = " << xhat << ", yhat = " << yhat << ", zhat = " << zhat << ". x.x = " << (xhat & xhat) << ", y.y = " << (yhat & yhat) <<", z.z = " << (zhat & zhat) << ", x.y = " << (xhat & yhat) << ", x.z = " << (xhat & zhat) << ", y.z = " << (yhat & zhat) << endl;
			pf1.append(pf1[0]);
        }
*/
        //Swapping pf1 points if pf0 and pf1 point in same general direction (because we want a quadrilateral ABCD where pf0 = AB and pf1 = CD)
        if ( ((B-A) & (D-C)) > 0 )
        {
			Info << "Swapping C and D" << endl;
            vector tmp = D;
            D = C;
            C = tmp;
        }

		Info << "A = " << A << ", B = " << B << ", C = " << C << ", D = " << D << endl;
        scalar Bx = mag(B-A);
        scalar Cx = (C-A) & xhat;
        scalar Cy = mag((C-A) & yhat);
        scalar Dx = (D-A) & xhat;
        scalar Dy = mag((D-A) & yhat);


//      area = ((Cx-Bx)*Dy-Dx*Cy)/6.0 + 0.25*Bx*(Dy+Cy);
        alpha = 0.5*((Cx-Bx)*Dy-Dx*Cy);
        beta = 0.5*Bx*(Dy+Cy);
		Info << "Bx = " << Bx << ", Cx = " << Cx << ", Cy = " << Cy << ", Dx = " << Dx << ", Dy = " << Dy << ", alpha = " << alpha << ", beta = " << beta << endl;
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
    scalarField& dVf,
	const scalar& dt
)
{
//    const scalarField& V = mesh_.V();
//    scalarField dV(alpha1_.size(),0.0);

//    const scalar tol = 1e-10;//SMALL;
    Info << "Bounding transport...\n" << endl;
	
	forAll(isSurfaceCell_,ci)
	{
		if (isSurfaceCell_[ci])
		{
			scalar waterGain = -netFlux(dVf,ci);
			scalar availableSpace = (1.0-alpha1_[ci])*(mesh_.V()[ci]);

			if (waterGain > availableSpace) //convert face air outfluxes to fractions of available air that is fluxed out
			{
				Info << "Cell " << ci << " with alpha1 = " << alpha1_[ci] << " has a waterGain = " << waterGain << " and availableSpace = " << availableSpace << "." << endl;
				scalar cellAirOutFlow(0.0);
				DynamicList<label> outFluxingFaces;
				getOutFluxFaces(ci, outFluxingFaces);				
				
				DynamicList<label> airOutFluxingFaces;
				DynamicList<scalar> faceAirOutFlows;
				forAll(outFluxingFaces,fi)
				{
					label fLabel = outFluxingFaces[fi];
					scalar faceAirOutFlow = mag(phi_[fLabel])*dt - mag(dVf[fLabel]);
					if (faceAirOutFlow > 0.0)
					{
						airOutFluxingFaces.append(fLabel);
						faceAirOutFlows.append(faceAirOutFlow);
						cellAirOutFlow += faceAirOutFlow;
					}
				}
				airOutFluxingFaces.shrink();
				faceAirOutFlows.shrink();
				Info << "Cell " << ci << ",\tairOutFluxingFaces = " << airOutFluxingFaces << endl; 
				if (cellAirOutFlow > 0.0)
				{
					Info << "cellAirOutFlow = " << cellAirOutFlow << endl;
					forAll(airOutFluxingFaces,fi)
					{
						label fLabel = airOutFluxingFaces[fi];
						scalar newAirOutFlow = availableSpace*(faceAirOutFlows[fi]/cellAirOutFlow);
						dVf[fLabel] = sign(dVf[fLabel])*(mag(phi_[fLabel])*dt-newAirOutFlow);
						Info << "Face " << fLabel << ": newAirOutFlow = " << newAirOutFlow << ", dVf = " << dVf[fLabel] << endl;
					}
					Info << "New waterGain = " << -netFlux(dVf,ci) << endl;
				}

			}
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
    const labelList fLabels = mesh_.cells()[cLabel];
    forAll(fLabels,fi)
    {
		const label fLabel = fLabels[fi];
        if (fLabel < mesh_.nInternalFaces())  //Error when this is removed - should be removed when boundary face treatment is implemented
        {
			const label owner = mesh_.owner()[fLabel];
            if (owner == cLabel)
            {
                dV += dVf[fLabel];
            }
            else
            {
                dV -= dVf[fLabel];
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
	Info << "dVf.size() = " << dVf.size() << ", mesh_.nFaces() = " << mesh_.nFaces() << endl;
    scalarField& dVfi = dVf;
    timeIntegratedFlux(dt, dVfi);
	boundAlpha(dVfi,dt);
	boundAlpha(dVfi,dt);
    alpha1_ -= fvc::surfaceIntegrate(dVf); //For each cell sum contributions from faces with pos sign for owner and neg sign for neighbour (as if it is a flux) and divide by cell volume
    alpha1_.correctBoundaryConditions();
//    alpha1_ = min(1.0,max(0.0,alpha1_));
//	scalar eps = 1e-10;
//    alpha1_ = alpha1_*pos(alpha1_-eps)*neg(alpha1_-(1.0-eps)) + pos(alpha1_-(1.0-eps));
}


// ************************************************************************* //
