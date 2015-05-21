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

void Foam::isoAdvection::timeIntegratedFlux
(
	const volScalarField& alpha,
	const surfaceScalarField& phi,
	const volVectorField& U,
	const scalar& dt,
	scalarField& dVf
)
{
	const fvMesh& mesh_ = alpha.mesh();
	volPointInterpolation vpi(mesh_);
	const scalarField alphap = vpi.interpolate(alpha);
    const cellList& cells = mesh_.cells();
	dVf = 0.0; //estimated total water volume transported across mesh faces during time interval dt

	isoCutter cutter(mesh_,alphap);
	
	forAll(cells,ci)
	{	
//		isoSubCell testCell(mesh_,ci);
		//Make list of all cell faces out of which fluid is flowing
		DynamicList<label> outFluxingFaces;
		getOutFluxFaces(phi, ci, outFluxingFaces);

		//Determine whether cell is fully above or fully below surface (if neither, it must be cut)
		scalar aMin, aMax;
		subSetExtrema(alphap,mesh_.cellPoints()[ci],aMin,aMax);

		//Set dVf on all outfluxing faces in accordance to the status of the cell (if allPointsAbove dVf should remain 0 on all outfluxing faces so do nothing)
//		if ( aMin > 0.5 )
		if ( aMin > 0.5 || aMax < 0.5 || alpha[ci] >= 1-SMALL || alpha[ci] <= SMALL )
		{
			forAll(outFluxingFaces,fi)
			{
				label lfi = outFluxingFaces[fi];
				if ( lfi < mesh_.nInternalFaces() )
				{
					dVf[lfi] = alpha[ci]*phi[lfi]*dt;
				}
			}
		}
		else if ( aMax > 0.5 ) //Then cell must be cut
		{
			Info << "\n------------ Cell " << ci << " ------------" << endl;
			Info << "outFluxingFaces: " << outFluxingFaces << endl;
			
			//Calculate isovalue f0 reproducing VOF value to given accuracy
			scalar f0, tol(1e-6);
			label maxIter(100);
			vector subCellCtr;
			cutter.vofCutCell(ci, alpha[ci], tol, maxIter, f0, subCellCtr);

			//Calculate isoFace0 centre xs0, normal ns0, and velocity U0 = U(xs0)
			vector x0(vector::zero), n0(vector::zero);
			cutter.isoFaceCentreAndArea(ci,f0,x0,n0); //Stupid to recalculate this here - should be provided by vofCutCell above

			//Make n0 point out of subCell
			Info << "x0: " << x0 << ", n0: " << n0 << ", subCellCtr: " << subCellCtr << endl;
			if (((x0 - subCellCtr) & n0) < -10*SMALL)
			{
				n0 *= (-1.0);
				Info << "Changing direction of n0 " << n0 << endl;
			}
			
			//Interpolate velocity to isoFace centre
			interpolationCellPoint<vector> UInterp(U);
			vector U0 = UInterp.interpolate(x0,ci);
			Info << "U0 = " << U0 << endl;
			scalar Un0 = (U0 & n0);
			Info << "Un0 = " << Un0 << endl;
			if (mag(n0) > SMALL)
			{
				Un0 /= mag(n0);
			}
			else
			{
				Info << "WARNING for cell " << ci << ": n0 has zero length. Probably because isoFace is a point. Skipping division of Un0 by mag(n0)." << endl;
				Info << "alpha = " << alpha[ci] << ", f0 = " << f0 << ", subCellCentre = " << subCellCtr << ", isoFaceCentre = " << x0 << ", isoFaceNormal = " << n0 << ", isoFaceVelocity = " << U0 << ", isoFaceNormalVelocity = " << Un0 << endl;
			}
			Info << "Initial values for time step:" << endl;
			Info << "f0 = " << f0 << ", subCellCentre = " << subCellCtr << ", isoFaceCentre = " << x0 << ", isoFaceNormal = " << n0 << ", isoFaceVelocity = " << U0 << ", isoFaceNormalVelocity = " << Un0 << endl;
			
			//Estimating time integrated water flux through each outfluxing face
			vector x00(x0), n00(n0); //Necessary when more than one outfluxing face
			scalar f00(f0);
			forAll(outFluxingFaces,fi)
			{
				label fLabel = outFluxingFaces[fi];
	
//------------------------Calculate initial alphaf on the face------------------------
				scalar t0 = 0.0; //Initial time
				scalar a0 = 0.0; //Initial alphaf
				x0 = x00; //Necessary when more than one outfluxing face
				n0 = n00; //Necessary when more than one outfluxing face
				f0 = f00;
				const vectorField& Sf = mesh_.Sf();
				pointField partSubFacePts;
				bool fullySubmerged = cutter.getSubFace(fLabel,f0,partSubFacePts);
				if (fullySubmerged)
				{
					a0 = 1;
				}
				else if (partSubFacePts.size() != 0)
				{
					vector fCtr, fArea;
					cutter.makeFaceCentreAndArea(partSubFacePts, fCtr, fArea);
					a0 = mag(fArea)/mag(Sf[fLabel]);
				}
				Info << "\nFACE " << fLabel << " with phi = " << phi[fLabel] << " and a0 = " << a0 << "." << endl;
				
				//Generating list of face vertices approached by isoFace
				scalarField fv; //List in ascending/descending order of unique face vertex alphap values larger/smaller than a0 for Un0 smaller/larger than zero (corresponding to backwards/forward motion of water surface)
				vector xFirst(vector::zero), xLast(vector::zero); //First and last vertex met by surface
				getVertexAlphas(alphap,mesh_,fLabel,f0,Un0,fv,xFirst,xLast);
				Info << "fv = " << fv << ", xFirst = " << xFirst << ", xLast = " << xLast << endl;
//----------------------------------------------------------------------
				
				if (fv.size() > 0) //((face is either completely submerged OR unsubmerged) AND Un is such that surface moves away from face) OR (face is cut and Un = 0).
				{
					label nv = 0;
					scalar a1(0.0); //VOF value on face corresponding to isoFace1
					scalar t1(0.0); //Next time a vertex is met
					while (t0 < dt && nv < fv.size())
					{
						Info << "Calculating flux integral contribution for sub time step #" << nv << " starting at t0 = " << t0 << " (dt = " << dt << ")." << endl;
						//Calculating new isoFace position and orientation
						scalar f1(fv[nv]);
						Info << "Calculating isoFace1 for the next met vertex alpha value, f1 = " << f1 << "." << endl;
						vector x1(vector::zero), n1(vector::zero); //IsoFace1 centre and area vectors
						
						if (nv == 0 && (a0 <= SMALL || a0 >= 1-SMALL)) //If face is initially completely full or empty so will it be at its first encounter with surface
						{
							a1 = a0;
							x1 = xFirst;
							n1 = n0;
							Info << "Since nv = " << nv << " and a0 = " << a0 << " we set a1 = a0 and x1 = xFirst = " << xFirst << endl;
						}
						else if (nv == fv.size()-1) //At the last point the face will either be completely filled or emptied - so isoFace = the vertex point
						{
							a1 = pos(f0-f1); //If f0 > f1, face is filling up with water - else it is being emptied
							x1 = xLast;
							n1 = n0; //since isoFace1 is a point its normal is undefined and we set it to previous value
							Info << "Since nv = " << nv << " = fv.size()-1 we set a1 = pos(f0-f1) = " << pos(f0-f1) << " and x1 = xLast = " << xLast << endl;
						} 
						else
						{
							cutter.isoFaceCentreAndArea(ci,f1,x1,n1);
	
							//Finding alphaf on face for new time
							fullySubmerged = cutter.getSubFace(fLabel,f0,partSubFacePts); //fLabel was previously fi - possible big bug!
							if (fullySubmerged)
							{
								a1 = 1;
							}
							else if (partSubFacePts.size() != 0)
							{
								vector fCtr, fArea;
								cutter.makeFaceCentreAndArea(partSubFacePts, fCtr, fArea);
								a1 = mag(fArea)/mag(Sf[fLabel]);
							}
						}
						Info << "New isoFace1 centre and area are x1 = " << x1 << " and n1 = " << n1 << " and results in a1 = " << a1 << " on face." << endl; 
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
//						Info << "x0: " << x0 << " x1: " << x1 << " n0: " << n0 << " U0: " << U0 << endl;
						t1 = min(t0 + dtn,dt); //Next time
						a1 = a0 + (t1-t0)/dtn*(a1-a0); //Linear interpolation of alphaf to time t1 from values at time t0 and t0+dtn - only changes a1 if t0+dtn > dt.
						f1 = f0 + (t1-t0)/dtn*(f1-f0);
						Info << "Estimated dtn = " << dtn << ", t1 = " << t1 << ", a1 = " << a1 << endl;
						dVf[fLabel] += phi[fLabel]*0.5*(a0 + a1)*(t1-t0); //Trapezoidal estimate
						Info << "dVf[fLabel] = " << dVf[fLabel] << " after adding " << phi[fLabel]*0.5*(a0 + a1)*(t1-t0) << " m3" << endl;
						scalar olddVf = phi[fLabel]*0.5*(a0 + a1)*(t1-t0);
						scalar newdVf = phi[fLabel]/mag(mesh_.Sf()[fLabel])*(t1-t0)*(a0*mag(mesh_.Sf()[fLabel]) + sign(Un0)*integratedArea(alphap,alpha,fLabel,f0,f1));
						if (mag(olddVf-newdVf) > 1e-10)
						{
							Info << "olddVf-newdVf: " << olddVf-newdVf << " m3" << endl;
						}
						t0 = t1;
						x0 = x1;
//						n0 = n1;
						a0 = a1;
						f0 = f1;
//						U0 = U1;
//						Un0 = Un1;
						nv++;
					}
					if (t1 < dt)
					{
						dVf[fLabel] += phi[fLabel]*a1*(dt-t1);
					}
				}
				else
				{
//					Info << "No elements in fVert. Face must be uncut and surface moving away from it." << endl;
					dVf[fLabel] += phi[fLabel]*a0*dt;
				}

			}
		}
		else
		{
			Info << "WARNING: Face not treated!!!" << endl;
		}
	}
}

void Foam::isoAdvection::getOutFluxFaces
(
	const surfaceScalarField& phi,
	const label& ci,
	DynamicList<label>& outFluxingFaces
)
{
	const fvMesh& mesh_ = phi.mesh();
    const cellList& cells = mesh_.cells();
	forAll(cells[ci],fi)
	{
		label fLabel = cells[ci][fi];
		if (mesh_.owner()[fLabel] == ci) 
		{
			if (phi[fLabel] > 10*SMALL) //HARDCODED LIMIT
			{
				outFluxingFaces.append(fLabel);
			}
		}
		else if ( phi[fLabel] < -10*SMALL ) //ci must be neighbour of fLabel
		{
			outFluxingFaces.append(fLabel);					
		}
	}
	outFluxingFaces.shrink();
}

void Foam::isoAdvection::getVertexAlphas
(
	const scalarField& alphap,
	const fvMesh& mesh_,
	const label& fLabel,
	const scalar& f0,
	const scalar& Un0,
	scalarField& fv,
	vector& xFirst,
	vector& xLast
)
{
	labelList pLabels = mesh_.faces()[fLabel];
//	Info << "pLabels: " << pLabels << endl;

	vectorField xVert; //Approached face vertex points
	scalarField fVert; //Corresponding alpha
	forAll(pLabels,pi)
	{
		scalar fVerti = alphap[pLabels[pi]];
//	Info << "fVerti for pi = " << pi << " equals " << fVerti << endl;
		if ( neg( Un0*(fVerti - f0) ) ) //Neg if Un0 and (fVerti-f0) have opposite signs i.e. if ( Un0 > 0 and (fVerti<f0) ) or ( Un0 < 0 and (fVerti>f0) )
		{
			xVert.append(mesh_.points()[pLabels[pi]]);
			fVert.append(fVerti);
		}
	}
//	Info << "fVert = " << fVert << ", xVert = " << xVert << endl;
	
	if (fVert.size() > 0) //((face is either completely submerged OR unsubmerged) AND Un is such that surface moves away from face) OR (face is cut and Un = 0).
	{
		scalar aMin(GREAT), aMax(-GREAT);
		label iaMin(-1), iaMax(-1);
		forAll(fVert,ai)
		{
			if (fVert[ai] < aMin)
			{
				aMin = fVert[ai];
				iaMin = ai;
			}
			if (fVert[ai] > aMax)
			{
				aMax = fVert[ai];
				iaMax = ai;
			}
		}
		xFirst = xVert[iaMax];
		xLast = xVert[iaMin];
		if (Un0 < 0)
		{
			Swap(xFirst,xLast);
		}
		sort(fVert);
		
//					Info << "Sorted fVert: " << fVert << endl;
		fv.append(fVert[0]);
		for (label n = 1; n < fVert.size(); n++)
		{
//			Info << "mag(fVert[n] - fVert[n-1]): " << mag(fVert[n] - fVert[n-1]) << endl;
			if (mag(fVert[n] - fVert[n-1]) > 10*SMALL) //SMALL is 1e-15 and I have experienced the difference being 1.2 e-15.
			{
				fv.append(fVert[n]);
			}
		}
//					Info << "Unique values of fVert in fv: " << fv << endl;
		if (Un0 > 0) //Water moves "forward", so nearest vertices are those with highest alphaf value
		{
			reverse(fv);
		}
	}
//	Info << "Reversed fv: " << fv << endl;	
}

Foam::scalar Foam::isoAdvection::integratedArea
(
	const scalarField& alphap,
	const volScalarField& alpha1,
	const label& fLabel,
	const scalar& f0,
	const scalar& f1
)
{
	const fvMesh& mesh = alpha1.mesh();
	isoCutter cutter(mesh,alphap);
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
		zhat = mesh.Sf()[fLabel]/mag(mesh.Sf()[fLabel]);
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
		
//		Info << "Bx = " << Bx << ", Cx = " << Cx << ", Cy = " << Cy << ", Dx = " << Dx << ", Dy = " << Dy << endl; 
		
		area = ((Cx-Bx)*Dy-Dx*Cy)/6.0 + 0.25*Bx*(Dy+Cy);
	}

	return area;
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
	const volScalarField& alpha,
	const surfaceScalarField& phi,
	scalarField& dVf
)
{
	const fvMesh& mesh_ = alpha.mesh();
	const scalarField& V = mesh_.V();
//	const cellList& cells = mesh_.cells(); //It is important for the efficiency of the code to verify that this call does not generate the cell faces for all cells each time it is called but only the first time!
//	const labelUList& owner = mesh_.owner();
	scalarField dV(alpha.size(),0.0);

	const scalar tol = 1e-8;//SMALL;
	Info << "Bounding transport (tolerance = " << tol << ")" << endl;
	
	DynamicList<label> overFilledCells, overEmptiedCells;
	forAll(alpha,ci)
	{
		dV[ci] = netFlux(dVf,mesh_,ci);
		bool cellIsOverFilled = dV[ci] < -tol && mag(dV[ci]) > V[ci]*(1.0-alpha[ci]) + tol;
		if ( cellIsOverFilled )
		{
//			Info << "Cell " << ci << " is overfilled: -dV = " << -dV[ci] << " < V*(1-alpha) = " << V[ci]*(1.0-alpha[ci]) << endl;
			overFilledCells.append(ci);			
		}
		else 
		{
			bool cellIsOverEmptied = dV[ci] > V[ci]*alpha[ci] + tol;
			if ( cellIsOverEmptied )
			{
//				Info << "Cell " << ci << " is overemptied: dV = " << dV[ci] << " > V*alpha = " << V[ci]*alpha[ci] << endl;
				overEmptiedCells.append(ci);
			}
		} 
	}
	overFilledCells.shrink();
	overEmptiedCells.shrink();

	//Correcting overfilling fluxes
	{
		label nOverFilledCells = overFilledCells.size();
		label ci = 0;
		while(ci < nOverFilledCells)
		{
			const label cLabel = overFilledCells[ci];
			DynamicList<label> outFluxingFaces;
			getOutFluxFaces(phi, cLabel, outFluxingFaces);
			scalar sumdVOut = 0.0;
			forAll(outFluxingFaces,fi) //If some outfluxing faces e.g. those fluxing pure air or water are to be left out, this is the place to do it.
			{
				const label fLabel = outFluxingFaces[fi];
				sumdVOut += mag(dVf[fLabel]);
			}
	//		scalar surplusWater = V[cLabel]*(alpha[cLabel]-1.0) - dV[cLabel];
			scalar surplusWater = mag(dV[cLabel])-V[cLabel]*(1.0-alpha[cLabel]);
			const scalar alphaNew = alpha[cLabel]-dV[cLabel]/V[cLabel];
			Info << "\nCell " << cLabel << " has alpha = " << alpha[cLabel] << ", dV = " << dV[cLabel] << " m3, surplusWater = " << surplusWater << " m3 which would give alphaNew = " << alphaNew << endl;
			if (surplusWater>tol)
			{
				Info << "sumdVOut = " << sumdVOut << endl;
				//			Info << "\nCell " << cLabel << " has surplusWater = " << surplusWater << " m3." << endl;
				forAll(outFluxingFaces,fi)
				{
					const label fLabel = outFluxingFaces[fi];
					scalar waterAdded = surplusWater*dVf[fLabel]/sumdVOut;
					dVf[fLabel] += waterAdded; //Whatever sign dVf has, it corresponds to out of cLabel. We must make it smaller, thus adding a little if it is negative and subtracting a little if it is positive
					//Checking if receiving cell becomes overfilled by the additional water it received
					label otherCellLabel;
					if (cLabel == mesh_.owner()[fLabel])
					{
						otherCellLabel = mesh_.neighbour()[fLabel];
					}
					else
					{
						otherCellLabel = mesh_.owner()[fLabel];	
					}
					scalar newdVOther = dV[otherCellLabel] - mag(waterAdded);
					bool otherCellIsOverfilled = newdVOther < -tol && mag(newdVOther) > V[otherCellLabel]*(1.0-alpha[otherCellLabel]) + tol;
					
					if ( otherCellIsOverfilled )
					{
						bool otherCellWasAlreadyOverfilled = dV[otherCellLabel] < -tol && (mag(dV[otherCellLabel]) > V[otherCellLabel]*(1.0-alpha[otherCellLabel]) + tol);
						if ( !otherCellWasAlreadyOverfilled )
						{
							Info << "Cell " << otherCellLabel << " was overfilled. Appending to overFilledCells list." << endl;
							overFilledCells.append(otherCellLabel);
							nOverFilledCells++;
						}
					}
					dV[otherCellLabel] = newdVOther; //Updating dV for neighbour cell - Do not move this in front of the above if loop
	//				dV[otherCellLabel] = netFlux(dVf,otherCellLabel);
					Info << "Moving " << waterAdded << "  m3 to its neighbour cell " << otherCellLabel << "." << endl;
				}
	//			dV[cLabel] = netFlux(dVf,cLabel);
				dV[cLabel] += surplusWater; //Updating dV for overfilled cell
			}
			ci++;
		}	
	}
////////////////////	
	//Correcting overemptying fluxes
	{
		label nOverEmptiedCells = overEmptiedCells.size();
		label ci = 0;
		while(ci < nOverEmptiedCells)
		{
			const label cLabel = overEmptiedCells[ci];
			DynamicList<label> outFluxingFaces;
			getOutFluxFaces(phi, cLabel, outFluxingFaces);
			scalar sumdVOut = 0.0;
			forAll(outFluxingFaces,fi) //If some outfluxing faces e.g. those fluxing pure air or water are to be left out, this is the place to do it.
			{
				const label fLabel = outFluxingFaces[fi];
				sumdVOut += mag(dVf[fLabel]);
			}
			scalar missingWater = dV[cLabel] - V[cLabel]*alpha[cLabel];
			const scalar alphaNew = alpha[cLabel]-dV[cLabel]/V[cLabel];
			Info << "\nCell " << cLabel << " has alpha = " << alpha[cLabel] << ", dV = " << dV[cLabel] << " m3, missingWater = " << missingWater << " m3 which would give alphaNew = " << alphaNew << endl;
			if (missingWater>tol)
			{
				Info << "sumdVOut = " << sumdVOut << endl;
				forAll(outFluxingFaces,fi)
				{
					const label fLabel = outFluxingFaces[fi];
					scalar waterSubtracted = missingWater*dVf[fLabel]/sumdVOut;
					dVf[fLabel] -= waterSubtracted; //Whatever sign dVf has, it corresponds to out of cLabel. We must make it smaller, thus adding a little if it is negative and subtracting a little if it is positive
					//Checking if receiving cell becomes overemptied by the decrease in water it receives
					label otherCellLabel;
					if (cLabel == mesh_.owner()[fLabel])
					{
						otherCellLabel = mesh_.neighbour()[fLabel];
					}
					else
					{
						otherCellLabel = mesh_.owner()[fLabel];	
					}
					scalar newdVOther = dV[otherCellLabel] + mag(waterSubtracted);
					bool otherCellIsOverEmptied = newdVOther > V[otherCellLabel]*alpha[otherCellLabel] + tol;
					if ( otherCellIsOverEmptied )
					{
						bool otherCellWasAlreadyOverEmptied = dV[otherCellLabel] > V[otherCellLabel]*alpha[otherCellLabel] + tol;
						if ( !otherCellWasAlreadyOverEmptied )
						{
							Info << "Cell " << otherCellLabel << " was overemptied. Appending to overEmptiedCells list." << endl;
							overEmptiedCells.append(otherCellLabel);
							nOverEmptiedCells++;
						}
					}
					dV[otherCellLabel] = newdVOther; //Updating dV for neighbour cell - Do not move this in front of the above if loop
	//				dV[otherCellLabel] = netFlux(dVf,otherCellLabel);
					Info << "Moving " << waterSubtracted << "  m3 from its neighbour cell " << otherCellLabel << "." << endl;
				}
				dV[cLabel] -= missingWater; //Updating dV for overemptied cell
	//			dV[cLabel] = netFlux(dVf,cLabel);
			}
			ci++;
		}	
	}
}

Foam::scalar Foam::isoAdvection::netFlux
(
	const scalarField& dVf,
	const fvMesh& mesh_,
	const label& cLabel
)
{
//	const fvMesh& mesh_ = dVf.mesh();
	scalar dV = 0.0;
	const labelList& fLabels = mesh_.cells()[cLabel];
	forAll(fLabels,fi)
	{
		if (fLabels[fi] < mesh_.nInternalFaces())  //Error when this is removed - should be removed when boundary face treatment is implemented
		{
			if (mesh_.owner()[fLabels[fi]] == cLabel)
			{
				dV += dVf[fLabels[fi]];
			}
			else
			{
				dV -= dVf[fLabels[fi]];				
			}
		}
	}	
	return dV;
}

void Foam::isoAdvection::advect
(
	volScalarField& alpha1,
	const surfaceScalarField& phi,
	const volVectorField& U,
	const scalar& dt
)
{
	surfaceScalarField dVf(0*phi); //Construct as copy - has same dimensions so will probably give problems. How to construct as zero's with phi's mesh etc e.g. surfaceScalarField dVf(phi.size(),0.0)?
	dVf.dimensions().reset(phi.mesh().V().dimensions());
	scalarField& dVfi = dVf;
	timeIntegratedFlux(alpha1, phi, U, dt, dVfi);
	boundAlpha(alpha1, phi, dVfi);
	alpha1 -= fvc::surfaceIntegrate(dVf); //For each cell sum contributions from faces with pos sign for owner and neg sign for neighbour (as if it is a flux) and divide by cell volume		
	alpha1.correctBoundaryConditions();
}


// ************************************************************************* //
