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

#ifdef ISODEBUG
#define isoDebug(x) x
#else
#define isoDebug(x) 
#endif 

Foam::isoAdvection::isoAdvection
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U,
    const dictionary& dict
)
:
    mesh_(alpha1.mesh()),
    alpha1_(alpha1),
    ap_(mesh_.nPoints(),0.0),
    phi_(phi),
    U_(U),
    isSurfaceCell_(mesh_.nCells(),false),
    nAlphaBounds_(dict.lookupOrDefault<label>("nAlphaBounds", 1)), 
    vof2IsoTol_(dict.lookupOrDefault<scalar>("vof2IsoTol", 1e-8)),
    surfCellTol_(dict.lookupOrDefault<scalar>("surfCellTol", 1e-8)),
    writeToLog_(dict.lookupOrDefault<bool>("writeToLog", true))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::isoAdvection::timeIntegratedFlux
(
    const scalar dt,
    surfaceScalarField& dVfa
)
{
    isoDebug(Info << "Enter timeIntegratedFlux" << endl;)

	//Estimated total water volume transported across mesh faces during time interval dt. The sign convenctino is like the flux phi, i.e. positive means out of owner cell.
	scalarField& dVf = dVfa.internalField();

    //Interpolating VOF field to mesh points
/*  
	//Testing 
	const labelListList& pCells = mesh_.pointCells();
    forAll(pCells,pi)
    {
        const labelList& pci = pCells[pi];
        ap_[pi] = 0.0;
        scalar vol(0.0);
        forAll(pci,ci)
        {
            const label cLabel = pci[ci];
            ap_[pi] += alpha1_[cLabel]*mesh_.V()[cLabel];
            vol += mesh_.V()[cLabel];
        }
        ap_[pi] /= vol;
    }
//	isoDebug(Info << "Using cell volume weighted cell-point interpolation" << endl;)
*/

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
            dVf[fi] = phi_[fi]*max(min(1,alpha1_[donorCell]),0)*dt;
        }
        else
        {
            isoDebug(Info << "Warning: No donor cell for face " << fi << endl;)
        }
    }

    //Isocutting estimate of water flux from surface cellss
    interpolationCellPoint<vector> UInterp(U_);
    forAll(surfaceCells,cellI)
    {
        const label ci = surfaceCells[cellI];
        isoDebug(Info << "\n------------ Cell " << ci << " with alpha1 = " << alpha1_[ci] << " and 1-alpha1 = " << 1.0-alpha1_[ci] << " ------------" << endl;)

        //Make list of all cell faces out of which fluid is flowing
        DynamicList<label> outFluxingFaces;
        getOutFluxFaces(ci, outFluxingFaces);
        isoDebug(Info << "outFluxingFaces: " << outFluxingFaces << "\n" << endl;)

        //Calculate isoFace0 centre xs0, normal ns0, and velocity U0 = U(xs0)
        scalar f0(0.0), Un0(0.0);
        vector x0(vector::zero), n0(vector::zero);
        calcIsoFace(ci,x0,n0,f0,Un0,UInterp); //This one really also should give us a0 on all faces since it is calculated anyway. Do this with a cutCell structure
        isoDebug(Info << "calcIsoFace gives initial surface: \nx0 = " << x0 << ", \nn0 = " << n0 << ", \nf0 = " << f0 << ", \nUn0 = " << Un0 << endl;)

        //Estimating time integrated water flux through each outfluxing face
        forAll(outFluxingFaces,fi)
        {
            const label fLabel = outFluxingFaces[fi];
            isoDebug(Info << "\nFace " << fLabel << " with outward normal nf = " << sign( (mesh_.Cf()[fLabel]-mesh_.C()[ci]) & mesh_.Sf()[fLabel] )*mesh_.Sf()[fLabel]/mag(mesh_.Sf()[fLabel]) << endl;)
            dVf[fLabel] = timeIntegratedFlux(fLabel,x0,n0,Un0,f0,dt);
        }
    }
}

void Foam::isoAdvection::findSurfaceCells
(
    DynamicList<label>& surfaceCells
)
{
    isoDebug(Info << "Enter findSurfaceCells" << endl;)

    isSurfaceCell_ = false;
    forAll(alpha1_,ci)
    {
//        scalar aMin(GREAT), aMax(-GREAT);
//        subSetExtrema(ap_,mesh_.cellPoints()[ci],aMin,aMax);
//        if ( (aMin < 0.5 && aMax > 0.5) || ( surfCellTol_ < alpha1_[ci] && alpha1_[ci] < 1.0-surfCellTol_ ) )
        if ( ( surfCellTol_ < alpha1_[ci] && alpha1_[ci] < 1.0-surfCellTol_ ) )
        {
            isSurfaceCell_[ci] = true;
            surfaceCells.append(ci);
        }
    }
    surfaceCells.shrink();
    Info << "\nnSurfaceCells = " << surfaceCells.size() << endl;
}

void Foam::isoAdvection::calcIsoFace
(
    const label ci,
    vector& x0,
    vector& n0,
    scalar& f0,
    scalar& Un0,
    const interpolationCellPoint<vector>& UInterp
)
{
    isoDebug(Info << "Enter calcIsoFace" << endl;)
    //Construction of isosurface calculator to get access to its functionality
    isoCutter cutter(mesh_,ap_);
    label maxIter(100);
    vector subCellCtr;
    cutter.vofCutCell(ci, alpha1_[ci], vof2IsoTol_, maxIter, f0, subCellCtr);
    cutter.isoFaceCentreAndArea(ci,f0,x0,n0); //Stupid to recalculate this here - should be provided by vofCutCell above
    isoDebug(Info << "f0 = " << f0 << ", x0 = " << x0 << ", n0 = " << n0 << endl;)

    if ( mag(n0) < 1.0e-8 ) //Cell almost full or empty so isoFace is undefined. Calculating normal by going a little into the cell
    {
        scalar aMin(GREAT), aMax(-GREAT);
		const labelList cellPts = mesh_.cellPoints()[ci];
        subSetExtrema(ap_,cellPts,aMin,aMax);
        if (mag(f0-aMin) < mag(f0-aMax)) //f0 almost equal to aMin, i.e. cell almost full
        {
            isoDebug(Info << "Cell is almost full with aMin = " << aMin << " and aMax = " << aMax << endl;)
            scalar f0Inside =  aMin + 1e-3*(aMax-aMin);
			vector xdum;
            cutter.isoFaceCentreAndArea(ci,f0Inside,xdum,n0);
            if (((x0 - mesh_.C()[ci]) & n0) < 0.0) //n0 should point from cell centre towards isoFace centre for an almost full cell
            {
                n0 *= (-1.0);
                isoDebug(Info << "Changing direction of n0 " << n0 << endl;)
            }
        }
        else if (mag(f0-aMin) > mag(f0-aMax))//f0 almost equal to aMax, i.e. cell almost empty
        {
            isoDebug(Info << "Cell is almost empty with aMin = " << aMin << " and aMax = " << aMax << endl;)
            scalar f0Inside =  aMax + 1e-3*(aMin-aMax);
			vector xdum;
            cutter.isoFaceCentreAndArea(ci,f0Inside,xdum,n0);
            if (((mesh_.C()[ci] - x0) & n0) < 0.0) //n0 should point from isoFace centre towards cell centre for an almost empty cell
            {
                n0 *= (-1.0);
                isoDebug(Info << "Changing direction of n0 " << n0 << endl;)
            }
        }
		else
		{
			Info << "Warning: aMax = aMin = " << aMax << " for cell " << ci << endl;
		}

    } //Make sure n0 points out of subCell
    else if (((x0 - subCellCtr) & n0) < 0.0) //n0 should point out of water surface, i.e. in the direction from subCell centre to isoFace centre
    {
        n0 *= (-1.0);
        isoDebug(Info << "Changing direction of n0 " << n0 << endl;)
    }

	if ( mag(n0) > 1.0e-8 )
	{
		isoDebug(Info << "Normalising n0: " << n0 << endl;)
		n0 /= mag(n0);
	}

    //Interpolate velocity to isoFace centre
    vector U0 = UInterp.interpolate(x0,ci);
    Un0 = (U0 & n0);
}

Foam::scalar Foam::isoAdvection::timeIntegratedFlux
(
    const label fLabel,
    const vector& x0,
    const vector& n0,
    const scalar Un0,
    const scalar f0,
    const scalar dt
)
{
    isoDebug(Info << "Enter timeIntegratedFlux for face " << fLabel << endl;)

    scalar dVf(0.0); //Volume flowing through face in time interval [0,dt] to be calculated below

    //Find sorted list of times where the isoFace will arrive at face points given initial position x0 and velocity Un0*n0
    const labelList& pLabels = mesh_.faces()[fLabel];
    const label nPoints = pLabels.size();
    pointField fPts(nPoints);
    forAll(fPts,pi)
    {
        fPts[pi] = mesh_.points()[pLabels[pi]];
    }
    scalarField pTimes(fPts.size());
	if ( mag(Un0) > 1e-12 )
	{
		pTimes = ((fPts - x0) & n0)/Un0; //Here we estimate time of arrival to the face points from their normal distance to the initial surface and the surface normal velocity
	}
	else //Un0 is almost zero
	{
        isoCutter cutter(mesh_,ap_);
		scalar alphaf;
		cutter.subFaceFraction(fLabel,f0,alphaf);
		dVf = phi_[fLabel]*dt*alphaf;
		return dVf;
	}

    scalarField sortedTimes(pTimes);
    sort(sortedTimes);
    isoDebug(Info << "sortedTimes = " << sortedTimes << endl;)

    //Dealing with case where face is not cut by surface during time interval [0,dt] because face was already passed by surface
    if ( sortedTimes[nPoints-1] <= 0.0 ) //All cuttings in the past
    {
        isoDebug(Info << "All cuttings in the past" << endl;)
        dVf = phi_[fLabel]*dt*pos(Un0); //If all face cuttings were in the past and cell is filling up (Un0>0) then face must be full during whole time interval
        return dVf;
    }

    //Dealing with case where face is not cut by surface during time interval [0,dt] because dt is too small for surface to reach closest face point
    if ( sortedTimes[0] >= dt ) //All cuttings in the future
    {
        isoDebug(Info << "All cuttings in the future" << endl;)
        dVf = phi_[fLabel]*dt*(1-pos(Un0)); //If all cuttings are in the future but non of them within [0,dt] then if cell is filling up (Un0 > 0) face must be empty during whole time interval
        return dVf;
    }

    //Cutting sortedTimes at 0 and dt
    DynamicList<scalar> t;
    t.append(0.0);
    forAll(sortedTimes,ti)
    {
        if ( 1e-3*dt < sortedTimes[ti] && sortedTimes[ti] < (1-1e-3)*dt )
        {
            t.append(sortedTimes[ti]);
        }
    }
    t.append(dt);
    isoDebug(Info << "Cutting sortedTimes at 0 and dt: t = " << t << endl;)

    bool faceUncutInFirstInterval(sortedTimes[0] > 0.0);
    bool faceUncutInLastInterval(sortedTimes[nPoints-1] < dt);

    //Dealing with cases where face is cut at least once during time interval [0,dt]
    label nt = 0; //Time interval counter
    scalar initialArea(0.0); //Submerged area at time
    if ( faceUncutInFirstInterval ) //Special treatment for first time interval if face is uncut during this
    {
        dVf = phi_[fLabel]*(t[nt+1]-t[nt])*(1.0-pos(Un0)); //If Un0 > 0 cell is filling up - hence if face is cut at a later time but not initially it must be initially empty
        initialArea = mag(mesh_.Sf()[fLabel])*(1.0-pos(Un0));
        isoDebug(Info << "faceUncutInFirstInterval, so special treatment for first time interval: [" << t[nt] << ", " << t[nt+1] << "] giving dVf = " << dVf << endl;)
        nt++;
    }
    else //calculate initialArea if face is initially cut
    {
        isoDebug(Info << "face is initially cut, so finding initial area, pTimes = " << pTimes << ", Un0 = " << Un0 << endl;)
        isoCutter cutter(mesh_,ap_);
        initialArea = cutter.getSubFaceFraction(fLabel, -sign(Un0)*pTimes, 0.0)*mag(mesh_.Sf()[fLabel]); //does not use ap_
//        scalar A(0.0), B(0.0);
//        isoDebug(Info << "Calculating quadCoeffs" << endl;)
//        quadAreaCoeffs(fLabel,pTimes,t[nt],t[nt+1],A,B);
    }

    isoDebug(Info << "InitialArea for next time step corresponds to face phase fraction a0 = " << initialArea/mag(mesh_.Sf()[fLabel]) << " where |Sf| = " << mag(mesh_.Sf()[fLabel]) << " was used." << endl;)
    while ( nt < t.size()-(1+faceUncutInLastInterval) ) //
    {
        scalar A(0.0), B(0.0);
        quadAreaCoeffs(fLabel,pTimes,t[nt],t[nt+1],A,B);
        scalar integratedQuadArea = sign(Un0)*(A/3.0 + 0.5*B); //Integration of area(t) = A*t^2+B*t from t = 0 to 1.
        scalar Unf = phi_[fLabel]/mag(mesh_.Sf()[fLabel]); //normal velocity on face
        dVf += Unf*(t[nt+1]-t[nt])*(initialArea + integratedQuadArea);
        initialArea += sign(Un0)*(A + B); //Adding quad area to submerged area
        isoDebug(Info << "Integrating area for " << nt+1 << "'th time interval: [" << t[nt] << ", " << t[nt+1] << "] giving dVf = " << dVf << " and a0 = " << initialArea/mag(mesh_.Sf()[fLabel]) << endl;)
//      Info << "face owner = " << mesh_.owner()[fLabel] << endl;
//      scalar newdVf = phi_[fLabel]/mag(mesh_.Sf()[fLabel])*(t1-t0)*(a0*mag(mesh_.Sf()[fLabel]) + sign(Un0)*integratedArea(fLabel,f0,f1));
        nt++;
    }

    if ( faceUncutInLastInterval ) //Special treatment for last time interval if face is uncut during this
    {
        isoDebug(Info << "faceUncutInLastInterval, so special treatment for last (" << nt+1 << "'th) time interval: [" << t[nt] << ", " << t[nt+1] << "]" << endl;)
        dVf += phi_[fLabel]*(t[nt+1]-t[nt])*pos(Un0); //If face is cut at some intermediate time but not at last time, then if Un0 > 0 (cell filling up) face must be filled at last time interval.
    }
    return dVf;
}


void Foam::isoAdvection::getOutFluxFaces
(
    const label ci,
    DynamicList<label>& outFluxingFaces
)
{
    const cellList& cells = mesh_.cells();
    const labelList& fLabels = cells[ci];
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
    const label fLabel,
    const label cLabel
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


void Foam::isoAdvection::quadAreaCoeffs
(
    const label fLabel,
    const scalarField& f,
    const scalar f0,
    const scalar f1,
    scalar& alpha,
    scalar& beta
)
{
    isoDebug(Info << "Enter quadAreaCoeffs" << endl;)
    isoCutter cutter(mesh_,ap_);
    DynamicList<point> pf0, pf1;
    cutter.getFaceCutPoints(fLabel,f,f0,pf0);
    cutter.getFaceCutPoints(fLabel,f,f1,pf1);

    label np0(pf0.size()), np1(pf1.size());
    isoDebug(Info << "Face " << fLabel << " was cut at " << pf0 << " by f0 = " << f0 << " and at " << pf1 << " by " << " f1 = " << f1 << endl;)

//  scalar area(0.0);
    alpha = 0.0;
    beta = 0.0;

    if ( np0 > 0 && np1 > 0)
    {
        //Defining quadrilateral vertices
        vector A(pf0[0]), C(pf1[0]), B(vector::zero), D(vector::zero);

        //Triangle cases
        if (np0 == 2)
        {
            B = pf0[1];
        }
        else if ( np0 == 1 )
        {
            B = A + 1e-4*(pf1[1]-pf1[0]);
        }
		else
		{
			Info << "Warning: Face " << fLabel << " was cut at " << pf0 << " by f0 = " << f0 << endl;			
		}

        if (np1 == 2)
        {
            D = pf1[1];
        }
        else if ( np1 == 1 )
        {
            D = C + 1e-4*(A-B);
        }
		else
		{
			Info << "Warning: Face " << fLabel << " was cut at " << pf1 << " by f1 = " << f1 << endl;			
		}

        //Defining local coordinates for area integral calculation
        vector zhat = mesh_.Sf()[fLabel]/mag(mesh_.Sf()[fLabel]);

        vector xhat = B-A;
        xhat -= (xhat & zhat)*zhat;
        xhat /= mag(xhat);

        vector yhat = zhat ^ xhat;
        yhat /= mag(yhat); //Should not be necessary

//        isoDebug(Info << "xhat = " << xhat << ", yhat = " << yhat << ", zhat = " << zhat << ". x.x = " << (xhat & xhat) << ", y.y = " << (yhat & yhat) <<", z.z = " << (zhat & zhat) << ", x.y = " << (xhat & yhat) << ", x.z = " << (xhat & zhat) << ", y.z = " << (yhat & zhat) << endl;)


        //Swapping pf1 points if pf0 and pf1 point in same general direction (because we want a quadrilateral ABCD where pf0 = AB and pf1 = CD)
        if ( ((B-A) & (D-C)) > 0 )
        {
////            Info << "Swapping C and D" << endl;
            vector tmp = D;
            D = C;
            C = tmp;
        }

////        Info << "A = " << A << ", B = " << B << ", C = " << C << ", D = " << D << endl;
        scalar Bx = mag(B-A);
        scalar Cx = (C-A) & xhat;
        scalar Cy = mag((C-A) & yhat);
        scalar Dx = (D-A) & xhat;
        scalar Dy = mag((D-A) & yhat);


//      area = ((Cx-Bx)*Dy-Dx*Cy)/6.0 + 0.25*Bx*(Dy+Cy);
        alpha = 0.5*((Cx-Bx)*Dy-Dx*Cy);
        beta = 0.5*Bx*(Dy+Cy);
////        Info << "Bx = " << Bx << ", Cx = " << Cx << ", Cy = " << Cy << ", Dx = " << Dx << ", Dy = " << Dy << ", alpha = " << alpha << ", beta = " << beta << endl;
        //area(t) = A*t^2+B*t
        //integratedArea = A/3+B/2
    }
	else
	{
		Info << "Face " << fLabel << " was cut at " << pf0 << " by f0 = " << f0 << " and at " << pf1 << " by " << " f1 = " << f1 << endl;
	}
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


void Foam::isoAdvection::limitFluxes
(
    surfaceScalarField& dVf,
    const scalar dt
)
{
	scalarField alphaNew = alpha1_ - fvc::surfaceIntegrate(dVf); 
	scalar aTol = 1.0e-12;
	scalar maxAlphaMinus1 = max(alphaNew-1);
	scalar minAlpha = min(alphaNew);

	for ( label n = 0; n < nAlphaBounds_; n++ )
    {
		Info << "Running bounding number " << n+1 << " of time " << mesh_.time().value() << endl;

		if ( maxAlphaMinus1 > aTol )
		{
			isoDebug(Info << "Bound from above... " << endl;)
			scalarField alpha1 = alpha1_.internalField();
//			scalarField dVfcorrected = dVf.internalField();
			scalarField dVfcorrected = dVf;
			DynamicList<label> correctedFaces;
			boundFromAbove(alpha1,dt,dVfcorrected,correctedFaces);
			forAll(correctedFaces,fi)
			{
				label fLabel = correctedFaces[fi];
				dVf[fLabel] = dVfcorrected[fLabel];
			}
		}

		if ( minAlpha < -aTol )
		{
			isoDebug(Info << "Bound from below... " << endl;)
			scalarField alpha2 = 1.0 - alpha1_.internalField();
			scalarField dVfcorrected = phi_.internalField()*dt - dVf;
//			dVfcorrected -= dVf; //phi_ and dVf have same sign and dVf is the portion of phi_*dt that is water. 
			//If phi_ > 0 then dVf > 0 and mag(phi_*dt-dVf) < mag(phi_*dt) as it should. 
			//If phi_ < 0 then dVf < 0 and mag(phi_*dt-dVf) < mag(phi_*dt) as it should. 
			DynamicList<label> correctedFaces;
			boundFromAbove(alpha2,dt,dVfcorrected,correctedFaces);
			forAll(correctedFaces,fi)
			{
				label fLabel = correctedFaces[fi];
				dVf[fLabel] = phi_[fLabel]*dt - dVfcorrected[fLabel];
			}
		}

		//Check if still unbounded
		alphaNew = alpha1_ - fvc::surfaceIntegrate(dVf); 
		maxAlphaMinus1 = max(alphaNew-1);
		minAlpha = min(alphaNew);
		scalar nUndershoots = sum(neg(alphaNew+aTol));
		scalar nOvershoots = sum(pos(alphaNew-1-aTol));
		Info << "After bounding number " << n+1 << " of time " << mesh_.time().value() << ":" << endl;
		Info << "nOvershoots = " << nOvershoots << " with max(alphaNew-1) = " << maxAlphaMinus1 << " and nUndershoots = " << nUndershoots << " with min(alphaNew) = " << minAlpha << endl;
	}
}


void Foam::isoAdvection::boundFromAbove
(
	const scalarField& alpha1,
	const scalar dt,
	scalarField& dVf,
	DynamicList<label>& correctedFaces
)
{
	correctedFaces.clear();
	
	scalar aTol = 1e-12;
//	scalar maxOvershoot = 0.0;
	forAll(alpha1,ci)
	{
		scalar Vi = mesh_.V()[ci];
		scalar alpha1New = alpha1[ci] - netFlux(dVf,ci)/Vi;
		scalar alphaOvershoot = alpha1New - 1.0;
/*		
		if ( alphaOvershoot > aTol )
		{
			scalar fluidToPassOn = alphaOvershoot*Vi;
			DynamicList<label> outFluxFaces;
			getOutFluxFaces(ci,outFluxFaces);
			scalar dVftot = 0.0;
			forAll(outFluxFaces,fi)
			{
				label fLabel = outFluxFaces[fi];
				dVftot += mag(phi_[fLabel]*dt);
			}
			forAll(outFluxFaces,fi)
			{
				label fLabel = outFluxFaces[fi];
				dVf[fLabel] += sign(phi_[fLabel])*fluidToPassOn*mag(phi_[fLabel]*dt)/dVftot;
				correctedFaces.append(fLabel);
			}
		}
*/		
		scalar fluidToPassOn = alphaOvershoot*Vi;
		label nFacesToPassFluidThrough = 1;
				
		//First try to pass surplus fluid on to neighbour cells that are not filled and to which dVf < phi*dt
		while ( alphaOvershoot > aTol && nFacesToPassFluidThrough > 0 )
		{
			isoDebug(Info << "\n\nBounding cell " << ci << " with alpha overshooting " << alphaOvershoot << endl;)
			//First find potential neighbour cells to pass surplus water to
			DynamicList<label> outFluxFaces;
			getOutFluxFaces(ci,outFluxFaces);
			DynamicList<label> facesToPassFluidThrough;
			DynamicList<scalar> dVfmax;
			scalar dVftot = 0.0;
			nFacesToPassFluidThrough = 0;
			
			isoDebug(Info << "OutFluxFaces: " << outFluxFaces << endl;)
			forAll(outFluxFaces,fi)
			{
				label fLabel = outFluxFaces[fi];
				scalar maxExtraFaceFluidTrans = mag(phi_[fLabel]*dt - dVf[fLabel]);
				//dVf has same sign as phi and so if phi>0 we have mag(phi_[fLabel]*dt) - mag(dVf[fLabel]) = phi_[fLabel]*dt - dVf[fLabel]
				//If phi<0 we have mag(phi_[fLabel]*dt) - mag(dVf[fLabel]) = -phi_[fLabel]*dt - (-dVf[fLabel]) > 0 since mag(dVf) < phi*dt
				isoDebug(Info << "outFluxFace " << fLabel << " has maxExtraFaceFluidTrans = " << maxExtraFaceFluidTrans << endl;)
				if ( maxExtraFaceFluidTrans/Vi > aTol )
//				if ( maxExtraFaceFluidTrans/Vi > aTol && mag(dVf[fLabel])/Vi > aTol ) //Last condition may be important because without this we will flux through uncut downwind faces
//				if ( true )
				{
					facesToPassFluidThrough.append(fLabel);
					dVfmax.append(maxExtraFaceFluidTrans);
					dVftot += mag(phi_[fLabel]*dt);
				}
			}

			isoDebug(Info << "\nfacesToPassFluidThrough: " << facesToPassFluidThrough << ", dVftot = " << dVftot << " m3 corresponding to dalpha = " << dVftot/Vi << endl;)
			forAll(facesToPassFluidThrough,fi)
			{
				label fLabel = facesToPassFluidThrough[fi];
				scalar fluidToPassThroughFace = fluidToPassOn*mag(phi_[fLabel]*dt)/dVftot;
				nFacesToPassFluidThrough += pos(dVfmax[fi] - fluidToPassThroughFace);
				fluidToPassThroughFace = min(fluidToPassThroughFace,dVfmax[fi]);
				label neiCell = otherCell(fLabel,ci);
				isoDebug(Info << "Passing " << fluidToPassThroughFace << " m3 through face " << fLabel << " to cell " << neiCell << endl;)
				isoDebug(Info << "This corresponds to lowering alpha of cell " << ci << " by " << fluidToPassThroughFace/Vi << endl;)
				isoDebug(Info << "dVf of face " << fLabel << " before correction: " << dVf[fLabel] << " m3 corresponding to da = " << dVf[fLabel]/Vi << endl;)
				dVf[fLabel] += sign(phi_[fLabel])*fluidToPassThroughFace;
				isoDebug(Info << "dVf of face " << fLabel << " after correction: " << dVf[fLabel] << " m3 corresponding to da = " << dVf[fLabel]/Vi << endl;)
				isoDebug(Info << "New neighbour cell alpha: " << alpha1[neiCell] - netFlux(dVf,neiCell)/mesh_.V()[neiCell] << endl;)
				correctedFaces.append(fLabel);
			}
			alpha1New = alpha1[ci] - netFlux(dVf,ci)/Vi;
			alphaOvershoot = alpha1New - 1.0;
			fluidToPassOn = alphaOvershoot*Vi;
			isoDebug(Info << "\nNew alpha for cell " << ci << ": " << alpha1New << endl;)
		}
/*		
		//Finally if any surplus fluid left pass it on to downwind neighbour cells
		if ( alphaOvershoot > aTol )
		{
			Info << "Could not pass on all surplus fluid from cell " << ci << " because all faces are already at their full flux capacity. Overshoot = " << alphaOvershoot << endl;
		}	
		
		if ( alphaOvershoot > maxOvershoot )
		{
			maxOvershoot = alphaOvershoot;
		}
*/
	}
	isoDebug(Info << "correctedFaces = " << correctedFaces << endl;)
}


Foam::scalar Foam::isoAdvection::netFlux
(
    const scalarField& dVf,
    const label cLabel
)
{
	//Net volume of water leaving cell cLabel. If dV is negative, cLabel, receives water from neighbours.
    scalar dV = 0.0; 
    const labelList& fLabels = mesh_.cells()[cLabel];
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
    const scalar dt
)
{
    isoDebug(Info << "Enter advect" << endl;)

	surfaceScalarField dVf
	(
		IOobject
		(
			"dVf",
			mesh_.time().timeName(),
			mesh_,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh_,
		dimensionedScalar("vol", dimVol, 0)
    );
    timeIntegratedFlux(dt, dVf);
	limitFluxes(dVf,dt);
	//For each cell sum contributions from faces with pos sign for owner and neg sign for neighbour (as if it is a flux) and divide by cell volume
    alpha1_ -= fvc::surfaceIntegrate(dVf); 
    alpha1_.correctBoundaryConditions();
}


// ************************************************************************* //