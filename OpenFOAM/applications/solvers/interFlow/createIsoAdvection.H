//IsoAdvector related stuff
isoAdvection advector(alpha1,phi,U);

//Read isoAdvector settings
const dictionary& isoAdvectorDict(mesh.solutionDict().subDict("isoAdvector"));
scalar snapAlphaTol(isoAdvectorDict.lookupOrDefault<scalar>("snapTol", 1e-12));
bool clipAlpha(isoAdvectorDict.lookupOrDefault<bool>("clip", true));
word interfaceMethod(isoAdvectorDict.lookupOrDefault<word>("interfaceMethod", "MULES"));

//Optional writing out of surface cell set
IOobject::writeOption writeOpt = IOobject::NO_WRITE;
const bool printSurfCells(isoAdvectorDict.lookupOrDefault<bool>("printSurfCells", false));
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

//Optional writing out of bounded cell set
writeOpt = IOobject::NO_WRITE;
const bool printBoundCells(isoAdvectorDict.lookupOrDefault<bool>("printBoundCells", false));

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

////Optional writing out of continuity errors
// writeOpt = IOobject::NO_WRITE;
// const bool writeContErrors(isoAdvectorDict.lookupOrDefault<bool>("writeContErrors", false));

// if (writeContErrors)
// {
    // writeOpt = IOobject::AUTO_WRITE;
// }
// volScalarField contErr
// (
   // IOobject
   // (
       // "contErr",
       // runTime.timeName(),
       // mesh,
       // IOobject::NO_READ,
       // writeOpt
   // ),
   // alpha1
// );

//Storage for unclipped alpha1 which will be written to file after clipped 
//version has been used in PIMPLE loop.
volScalarField alpha1Unbounded = alpha1;

//Preparing for usage with prescribed and optionally periodic or reversed U
bool prescribedU(isoAdvectorDict.lookupOrDefault<bool>("prescribedU", false));
const scalar period = isoAdvectorDict.lookupOrDefault<scalar>("period", 0.0);
scalar reverseTime = isoAdvectorDict.lookupOrDefault<scalar>("reverseTime", 0.0);

//Creating copies of initial U and phi for periodic flows
//Rewrite to  avoid these extra fields when predcribedU = false
volVectorField U0 = U;
surfaceScalarField phi0 = phi;


//Timing advection step for log file
scalar advectionTime = 0;

//Initial Volume of Fluid (VOF)
const scalar V0(gSum(mesh.V().field()*alpha1.internalField()) + SMALL);
