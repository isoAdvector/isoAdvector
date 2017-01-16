//IsoAdvector related stuff
isoAdvection advector(alpha1,phi,U);

//Initial Volume of Fluid (VOF)
const scalar V0 = gSum(mesh.V().field()*alpha1.internalField()) + SMALL;

//Read isoAdvector settings
const dictionary& isoAdvectorDict = mesh.solutionDict().subDict("isoAdvector");
scalar clipAlphaTol = isoAdvectorDict.lookupOrDefault<scalar>("clipAlphaTol", 0);
bool snapAlpha = isoAdvectorDict.lookupOrDefault<bool>("snapAlpha", true);
bool useIsoAdvector = isoAdvectorDict.lookupOrDefault<bool>("useIsoAdvector", false);

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
//volScalarField contErr
//(
//    IOobject
//    (
//        "contErr",
//        runTime.timeName(),
//        mesh,
//        IOobject::NO_READ,
//        IOobject::AUTO_WRITE
//    ),
//    alpha1
//);

//Storage for unclipped alpha1 which will be written to file after clipped 
//version has been used in PIMPLE loop.
volScalarField alpha1Unbounded = alpha1;

scalar advectionTime = 0;

