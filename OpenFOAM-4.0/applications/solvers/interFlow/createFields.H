Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

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


Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();


// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);


// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

// MULES flux from previous time-step
surfaceScalarField alphaPhi
(
    IOobject
    (
        "alphaPhi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    phi*fvc::interpolate(alpha1)
);

// MULES Correction
//tmp<surfaceScalarField> talphaPhiCorr0;

#include "createMRF.H"

//IsoAdvector related stuff
isoAdvection advector(alpha1,phi,U);

//Initial Volume of Fluid (VOF)
const scalar V0 = gSum(mesh.V().field()*alpha1.internalField());

//Read isoAdvector settings
const dictionary& isoAdvectorDict = mesh.solutionDict().subDict("isoAdvector");
const scalar clipAlphaTol = isoAdvectorDict.lookupOrDefault<scalar>("clipAlphaTol", 0);
const bool snapAlpha = isoAdvectorDict.lookupOrDefault<bool>("snapAlpha", true);

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

//Solver writes out continuity error field
volScalarField contErr
(
    IOobject
    (
        "contErr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha1
);

//Storage for unclipped alpha1 which will be written to file after clipped 
//version has been used in PIMPLE loop.
volScalarField alpha1Org = alpha1;