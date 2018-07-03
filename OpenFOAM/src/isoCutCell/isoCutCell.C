/*---------------------------------------------------------------------------*\
              Original work | Copyright (C) 2016-2017 DHI
              Modified work | Copyright (C) 2016-2017 OpenCFD Ltd.
              Modified work | Copyright (C) 2017-2018 Johan Roenby
-------------------------------------------------------------------------------

License
    This file is part of isoAdvector which is an extension to OpenFOAM.

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
#include "scalarMatrices.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Debugging * * * * * * * * * * * * * //

#ifndef DebugPout
//Taken from OpenFOAM-plus/src/OpenFOAM/db/error/messageStream.H to make code
//compile with older OF versions.
#define DebugPout                                                              \
    if (debug) Pout
#endif

#ifndef InfoInFunction
#define InfoInFunction InfoIn(__func__)
#endif

#ifndef DebugInFunction
#define DebugInFunction                                                        \
    if (debug) InfoInFunction
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::isoCutCell::debug = 0;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isoCutCell::isoCutCell(const fvMesh& mesh, scalarField& f)
:
    mesh_(mesh),
    cellI_(-1),
    f_(f),
    isoValue_(0),
    isoCutFace_(isoCutFace(mesh_, f_)),
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
    cellStatus_(-1),
    subCellCentreAndVolumeCalculated_(false),
    isoFaceCentreAndAreaCalculated_(false)
{
    clearStorage();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::isoCutCell::calcSubCellCentreAndVolume()
{
    if (cellStatus_ == 0)
    {
        subCellCentre_ = vector::zero;
        subCellVolume_ = 0.0;

        // Estimate the approximate cell centre as the average of face centres
        label nCellFaces(1 + isoCutFaceCentres_.size() + fullySubFaces_.size());
        vector cEst = isoFaceCentre_ + sum(isoCutFaceCentres_);
        forAll(fullySubFaces_, facei)
        {
            cEst += mesh_.faceCentres()[fullySubFaces_[facei]];
        }
        cEst /= scalar(nCellFaces);


        // Contribution to subcell centre and volume from isoface
        const scalar pyr3Vol0 =
            max(mag(isoFaceArea_ & (isoFaceCentre_ - cEst)), VSMALL);

        // Calculate face-pyramid centre
        const vector pc0 = 0.75*isoFaceCentre_ + 0.25*cEst;

        // Accumulate volume-weighted face-pyramid centre
        subCellCentre_ += pyr3Vol0*pc0;

        // Accumulate face-pyramid volume
        subCellVolume_ += pyr3Vol0;

        // Contribution to subcell centre and volume from cut faces
        forAll(isoCutFaceCentres_, facei)
        {
            // Calculate 3*face-pyramid volume
            scalar pyr3Vol =
                max
                (
                    mag
                    (
                        isoCutFaceAreas_[facei]
                      & (isoCutFaceCentres_[facei] - cEst)
                    ),
                    VSMALL
                );

            // Calculate face-pyramid centre
            vector pc = 0.75*isoCutFaceCentres_[facei] + 0.25*cEst;

            // Accumulate volume-weighted face-pyramid centre
            subCellCentre_ += pyr3Vol*pc;

            // Accumulate face-pyramid volume
            subCellVolume_ += pyr3Vol;
        }

        // Contribution to subcell centre and volume from fully submerged faces
        forAll(fullySubFaces_, i)
        {
            const label facei = fullySubFaces_[i];
            const point& fCentre = mesh_.faceCentres()[facei];
            const vector& fArea = mesh_.faceAreas()[facei];

            // Calculate 3*face-pyramid volume
            scalar pyr3Vol = max(mag(fArea & (fCentre - cEst)), VSMALL);

            // Calculate face-pyramid centre
            vector pc = 0.75*fCentre + 0.25*cEst;

            // Accumulate volume-weighted face-pyramid centre
            subCellCentre_ += pyr3Vol*pc;

            // Accumulate face-pyramid volume
            subCellVolume_ += pyr3Vol;
        }

        subCellCentre_ /= subCellVolume_;
        subCellVolume_ /= scalar(3);
        VOF_ = subCellVolume_/mesh_.cellVolumes()[cellI_];

        subCellCentreAndVolumeCalculated_ = true;

        if (debug)
        {
            vector sumSf = isoFaceArea_;
            scalar sumMagSf = mag(isoFaceArea_);
            forAll(isoCutFaceCentres_, facei)
            {
                sumSf += isoCutFaceAreas_[facei];
                sumMagSf += mag(isoCutFaceAreas_[facei]);
            }
            forAll(fullySubFaces_, facei)
            {
                sumSf += mesh_.faceAreas()[fullySubFaces_[facei]];
                sumMagSf += mag(isoCutFaceAreas_[facei]);
            }
            if (mag(sumSf) > 1e-10)
            {
                Pout<< "Warninig: mag(sumSf)/magSumSf = "
                    << mag(sumSf)/sumMagSf << " for surface cell"
                    << cellI_ << endl;
            }
        }
    }
    else if (cellStatus_ == 1)
    {
        // Cell fully above isosurface
        subCellCentre_ = vector::zero;
        subCellVolume_ = 0;
        VOF_ = 0;
    }
    else if (cellStatus_ == -1)
    {
        // Cell fully below isosurface
        subCellCentre_ = mesh_.cellCentres()[cellI_];
        subCellVolume_ = mesh_.cellVolumes()[cellI_];
        VOF_ = 1;
    }
}


void Foam::isoCutCell::calcIsoFaceCentreAndArea()
{
    // Initial guess of face centre from edge points
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

    if (nEdgePoints > 0)
    {
        fCentre /= nEdgePoints;
    }
    else
    {
        DebugPout << "Warning: nEdgePoints = 0 for cell " << cellI_ << endl;
    }

    vector sumN = vector::zero;
    scalar sumA = 0;
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

            // Edges may have different orientation
            sumN += Foam::sign(n & sumN)*n;
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
        isoFaceCentre_ = sumAc/sumA/scalar(3);
        isoFaceArea_ = 0.5*sumN;
    }


    // Check isoFaceArea_ direction and change if not pointing out of subcell
    if ((isoFaceArea_ & (isoFaceCentre_ - subCellCentre())) < 0)
    {
        isoFaceArea_ *= (-1);
    }

    isoFaceCentreAndAreaCalculated_ = true;
}


void Foam::isoCutCell::calcIsoFacePointsFromEdges()
{
    DebugPout
        << "Enter calcIsoFacePointsFromEdges() with isoFaceArea_ = "
        << isoFaceArea_ << " and isoFaceCentre_ = " << isoFaceCentre_
        << " and isoFaceEdges_ = " << isoFaceEdges_ << endl;

    // Defining local coordinates with zhat along isoface normal and xhat from
    // isoface centre to first point in isoFaceEdges_
    const vector zhat = isoFaceArea_/mag(isoFaceArea_);
    vector xhat = isoFaceEdges_[0][0] - isoFaceCentre_;
    xhat = (xhat - (xhat & zhat)*zhat);
    xhat /= mag(xhat);
    vector yhat = zhat^xhat;
    yhat /= mag(yhat);

    DebugPout << "Calculated local coordinates" << endl;

    // Calculating isoface point angles in local coordinates
    DynamicList<point> unsortedIsoFacePoints(3*isoFaceEdges_.size());
    DynamicList<scalar> unsortedIsoFacePointAngles(3*isoFaceEdges_.size());
    forAll(isoFaceEdges_, ei)
    {
        const DynamicList<point>& edgePoints = isoFaceEdges_[ei];
        forAll(edgePoints, pi)
        {
            const point& p = edgePoints[pi];
            unsortedIsoFacePoints.append(p);
            unsortedIsoFacePointAngles.append
            (
                Foam::atan2
                (
                    ((p - isoFaceCentre_) & yhat),
                    ((p - isoFaceCentre_) & xhat)
                )
            );
        }
    }

    DebugPout<< "Calculated isoFace point angles" << endl;

    // Sorting isoface points by angle and inserting into isoFacePoints_
    labelList order(unsortedIsoFacePointAngles.size());
    Foam::sortedOrder(unsortedIsoFacePointAngles, order);
    isoFacePoints_.append(unsortedIsoFacePoints[order[0]]);
    for (label pi = 1; pi < order.size(); pi++)
    {
        if
        (
            // Note: Hardcoded tolerance
            mag
            (
                unsortedIsoFacePointAngles[order[pi]]
               -unsortedIsoFacePointAngles[order[pi-1]]
            ) > 1e-8
        )
        {
            isoFacePoints_.append(unsortedIsoFacePoints[order[pi]]);
        }
    }

    DebugPout<< "Sorted isoface points by angle" << endl;
}


Foam::label Foam::isoCutCell::calcSubCell
(
    const label celli,
    const scalar isoValue
)
{
    // Populate isoCutFaces_, isoCutFacePoints_, fullySubFaces_, isoFaceCentre_
    // and isoFaceArea_.

    clearStorage();
    cellI_ = celli;
    isoValue_ = isoValue;
    const cell& c = mesh_.cells()[celli];

    forAll(c, fi)
    {
        const label facei = c[fi];

        const label faceStatus = isoCutFace_.calcSubFace(facei, isoValue_);

        if (faceStatus == 0)
        {
            // Face is cut
            isoCutFacePoints_.append(isoCutFace_.subFacePoints());
            isoCutFaceCentres_.append(isoCutFace_.subFaceCentre());
            isoCutFaceAreas_.append(isoCutFace_.subFaceArea());
            isoFaceEdges_.append(isoCutFace_.surfacePoints());
        }
        else if (faceStatus == -1)
        {
            // Face fully below
            fullySubFaces_.append(facei);
        }
    }

    if (isoCutFacePoints_.size())
    {
        // Cell cut at least at one face
        cellStatus_ = 0;
        calcIsoFaceCentreAndArea();

        // In the rare but occuring cases where a cell is only touched at a
        // point or a line the isoFaceArea_ will have zero length and here the
        // cell should be treated as either completely empty or full.
        if (mag(isoFaceArea_) < 10*SMALL)
        {
            if (fullySubFaces_.empty())
            {
                // Cell fully above isosurface
                cellStatus_ = 1;
            }
            else
            {
                // Cell fully below isosurface
                cellStatus_ = -1;
            }
        }
    }
    else if (fullySubFaces_.empty())
    {
        // Cell fully above isosurface
        cellStatus_ = 1;
    }
    else
    {
        // Cell fully below isosurface
        cellStatus_ = -1;
    }

    return cellStatus_;
}


const Foam::point& Foam::isoCutCell::subCellCentre()
{
    if (!subCellCentreAndVolumeCalculated_)
    {
        calcSubCellCentreAndVolume();
    }

    return subCellCentre_;
}


Foam::scalar Foam::isoCutCell::subCellVolume()
{
    if (!subCellCentreAndVolumeCalculated_)
    {
        calcSubCellCentreAndVolume();
    }

    return subCellVolume_;
}


const Foam::DynamicList<Foam::point>& Foam::isoCutCell::isoFacePoints()
{
    if (cellStatus_ == 0 && isoFacePoints_.size() == 0)
    {
        calcIsoFacePointsFromEdges();
    }

    return isoFacePoints_;
}


const Foam::point& Foam::isoCutCell::isoFaceCentre()
{
    if (!isoFaceCentreAndAreaCalculated_)
    {
        calcIsoFaceCentreAndArea();
    }

    return isoFaceCentre_;
}


const Foam::vector& Foam::isoCutCell::isoFaceArea()
{
    if (!isoFaceCentreAndAreaCalculated_)
    {
        calcIsoFaceCentreAndArea();
    }

    return isoFaceArea_;
}


Foam::scalar Foam::isoCutCell::volumeOfFluid()
{
    if (!subCellCentreAndVolumeCalculated_)
    {
        calcSubCellCentreAndVolume();
    }

    return VOF_;
}


Foam::scalar Foam::isoCutCell::isoValue() const
{
    return isoValue_;
}


void Foam::isoCutCell::clearStorage()
{
    cellI_ = -1;
    isoValue_ = 0;
    isoCutFace_.clearStorage();
    isoCutFaces_.clear();
    isoCutFacePoints_.clear();
    isoCutFaceCentres_.clear();
    isoCutFaceAreas_.clear();
    isoFaceEdges_.clear();
    isoFacePoints_.clear();
    isoFaceCentre_ = vector::zero;
    isoFaceArea_ = vector::zero;
    subCellCentre_ = vector::zero;
    subCellVolume_ = -10;
    VOF_ = -10;
    fullySubFaces_.clear();
    cellStatus_ = -1;
    subCellCentreAndVolumeCalculated_ = false;
    isoFaceCentreAndAreaCalculated_ = false;
}


Foam::label Foam::isoCutCell::vofCutCell
(
    const label celli,
    const scalar alpha1,
    const scalar tol,
    const label maxIter
)
{
    DebugInFunction
        << "vofCutCell for cell " << celli << " with alpha1 = "
        << alpha1 << " ------" << endl;

    // Finding cell vertex extremum values
    const labelList& pLabels = mesh_.cellPoints(celli);
    scalarField fvert(pLabels.size());
    forAll(pLabels, pi)
    {
        fvert[pi] = f_[pLabels[pi]];
    }
    labelList order(fvert.size());
    sortedOrder(fvert, order);
    scalar f1 = fvert[order.first()];
    scalar f2 = fvert[order.last()];

    DebugPout << "fvert = " << fvert << ", and order = " << order << endl;

    // Handling special case where method is handed an almost full/empty cell
    if (alpha1 < tol)
    {
        return calcSubCell(celli, f2);
    }
    else if (1 - alpha1 < tol)
    {
        return calcSubCell(celli, f1);
    }

    // Finding the two vertices inbetween which the isovalue giving alpha1 lies
    label L1 = 0;
    label L2 = fvert.size() - 1;
    scalar a1 = 1;
    scalar a2 = 0;
    scalar L3, f3, a3;

    while (L2 - L1 > 1)
    {
        L3 = round(0.5*(L1 + L2));
        f3 = fvert[order[L3]];
        calcSubCell(celli, f3);
        a3 = volumeOfFluid();
        if (a3 > alpha1)
        {
            L1 = L3; f1 = f3; a1 = a3;
        }
        else if (a3 < alpha1)
        {
            L2 = L3; f2 = f3; a2 = a3;
        }
    }

    if (mag(f1 - f2) < 10*SMALL)
    {
        DebugPout<< "Warning: mag(f1 - f2) < 10*SMALL" << endl;
        return calcSubCell(celli, f1);
    }

    if (mag(a1 - a2) < tol)
    {
        DebugPout<< "Warning: mag(a1 - a2) < tol for cell " << celli << endl;
        return calcSubCell(celli, 0.5*(f1 + f2));
    }

    // Now we know that a(f) = alpha1 is to be found on the f interval
    // [f1, f2], i.e. alpha1 will be in the interval [a2,a1]
    DebugPout
        << "L1 = " << L1 << ", f1 = " << f1 << ", a1 = " << a1 << nl
        << "L2 = " << L2 << ", f2 = " << f2  << ", a2 = " << a2 << endl;


    // Finding coefficients in 3 deg polynomial alpha(f) from 4 solutions

    // Finding 2 additional points on 3 deg polynomial
    f3 = f1 + (f2 - f1)/scalar(3);
    calcSubCell(celli, f3);
    a3 = volumeOfFluid();

    scalar f4 = f1 + (f2 - f1)*scalar(2)/scalar(3);
    calcSubCell(celli, f4);
    scalar a4 = volumeOfFluid();

    // Building and solving Vandermonde matrix equation
    scalarField a(4), f(4), C(4);
    {
        a[0] = a1, f[0] = 0;
        a[1] = a3, f[1] = (f3 - f1)/(f2 - f1);
        a[2] = a4, f[2] = (f4 - f1)/(f2 - f1);
        a[3] = a2, f[3] = 1;
        scalarSquareMatrix M(4);
        forAll(f, i)
        {
            forAll(f, j)
            {
                M[i][j] = pow(f[i], 3 - j);
            }
        }

        // C holds the 4 polynomial coefficients
        C = a;
        LUsolve(M, C);
    }

    // Finding root with Newton method
    f3 = f[1]; a3 = a[1];
    label nIter = 0;
    scalar res = mag(a3 - alpha1);
    while (res > tol && nIter < 10*maxIter)
    {
        f3 -=
            (C[0]*pow3(f3) + C[1]*sqr(f3) + C[2]*f3 + C[3] - alpha1)
           /(3*C[0]*sqr(f3) + 2*C[1]*f3 + C[2]);
        a3 = C[0]*pow3(f3) + C[1]*sqr(f3) + C[2]*f3 + C[3];
        res = mag(a3 - alpha1);
        nIter++;
    }
    // Scaling back to original range
    f3 = f3*(f2 - f1) + f1;

    // Check result
    calcSubCell(celli, f3);
    const scalar VOF = volumeOfFluid();
    res = mag(VOF - alpha1);

    if (res > tol)
    {
        DebugPout
            << "Newton obtained f3 = " << f3 << " and a3 = " << a3
            << " with mag(a3-alpha1) = " << mag(a3-alpha1)
            << " but calcSubCell(celli,f3) gives VOF  = " << VOF << nl
            << "M(f)*C = a with " << nl
            << "f_scaled = " << f << nl
            << "f = " << f*(f2 - f1) + f1 << nl
            << "a = " << a << nl
            << "C = " << C << endl;
    }
    else
    {
        DebugPout<< "Newton did the job" << endl;
        return cellStatus_;
    }

    // If tolerance not met use the secant method  with f3 as a hopefully very
    // good initial guess to crank res the last piece down below tol
    // Note: This is expensive because subcell is recalculated every iteration
    scalar x2 = f3;
    scalar g2 = VOF - alpha1;
    scalar x1 = max(1e-3*(f2 - f1), 100*SMALL);
    x1 = min(max(x1, f1), f2);
    calcSubCell(celli, x1);
    scalar g1 = volumeOfFluid() - alpha1;

    nIter = 0;
    scalar g0(0), x0(0);
    while (res > tol && nIter < maxIter && g1 != g2)
    {
        x0 = (x2*g1 - x1*g2)/(g1 - g2);
        calcSubCell(celli, x0);
        g0 = volumeOfFluid() - alpha1;
        res = mag(g0);
        x2 = x1; g2 = g1;
        x1 = x0; g1 = g0;
        nIter++;
    }

    if (debug)
    {
        if (res < tol)
        {
            Pout<< "Bisection finished the job in " << nIter << " iterations."
                << endl;
        }
        else
        {
            Pout<< "Warning: Bisection not converged " << endl;
            Pout<< "Leaving vofCutCell with f3 = " << f3 << " giving a3 = "
                << a3 << " so alpha1 - a3 = " << alpha1 - a3 << endl;
        }
    }

    return cellStatus_;
}


void Foam::isoCutCell::volumeOfFluid
(
    volScalarField& alpha1,
    const scalar f0
)
{
    // Setting internal field
    scalarField& alphaIn = alpha1;
    forAll(alphaIn, celli)
    {
        const label cellStatus = calcSubCell(celli, f0);
        if (cellStatus != 1)
        {
            // If cell not entirely above isosurface
            alphaIn[celli] = volumeOfFluid();
        }
    }

    // Setting boundary alpha1 values
    forAll(mesh_.boundary(), patchi)
    {
        if (mesh_.boundary()[patchi].size() > 0)
        {
            const label start = mesh_.boundary()[patchi].patch().start();
            scalarField& alphap = alpha1.boundaryFieldRef()[patchi];
            const scalarField& magSfp = mesh_.magSf().boundaryField()[patchi];

            forAll(alphap, patchFacei)
            {
                const label facei = patchFacei + start;
                const label faceStatus = isoCutFace_.calcSubFace(facei, f0);

                if (faceStatus != 1)
                {
                    // Face not entirely above isosurface
                    alphap[patchFacei] =
                        mag(isoCutFace_.subFaceArea())/magSfp[patchFacei];
                }
            }
        }
    }
}


// ************************************************************************* //
