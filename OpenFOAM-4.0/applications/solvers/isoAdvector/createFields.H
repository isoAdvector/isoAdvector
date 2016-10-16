Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

Info<< "Reading field alpha1\n" << endl;
volScalarField alpha1
(
    IOobject
    (
        "alpha1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//IsoAdvector related stuff
isoAdvection advector(alpha1,phi,U);

//Initial Volume of Fluid (VOF)
const scalar V0 = gSum(mesh.V().field()*alpha1.internalField());

//Read isoAdvector settings
const dictionary& isoAdvectorDict = mesh.solutionDict().subDict("isoAdvector");
const scalar clipAlphaTol = isoAdvectorDict.lookupOrDefault<scalar>("clipAlphaTol", 0);
const bool snapAlpha = isoAdvectorDict.lookupOrDefault<bool>("snapAlpha", true);
const scalar period = isoAdvectorDict.lookupOrDefault<scalar>("period", 0.0);
scalar reverseTime = isoAdvectorDict.lookupOrDefault<scalar>("reverseTime", 0.0);

//Copies of initial U and phi for use when flow is periodic
surfaceScalarField phi0 = phi;
volVectorField U0 = U;

//Optional writing out of surface cell set and bounded cell sets
IOobject::writeOption writeOpt = IOobject::NO_WRITE;
const bool printSurfCells
    = isoAdvectorDict.lookupOrDefault<bool>("printSurfCells", false);
if (printSurfCells)
{
    writeOpt = IOobject::AUTO_WRITE;
}

cellSet surfCells
(
    IOobject
    (
        "surfCells",
        runTime.timeName(),
        runTime,
        IOobject::NO_READ,
        writeOpt
    )
);

writeOpt = IOobject::NO_WRITE;
const bool printBoundCells
    = isoAdvectorDict.lookupOrDefault<bool>("printBoundCells", false);

if (printBoundCells)
{
    writeOpt = IOobject::AUTO_WRITE;
}
cellSet boundCells
(
    IOobject
    (
        "boundCells",
        runTime.timeName(),
        runTime,
        IOobject::NO_READ,
        writeOpt
    )
);
