const dictionary& isoAdvectorControl = mesh.solutionDict().subDict("isoAdvector");


const scalar vof2IsoTol = isoAdvectorControl.lookupOrDefault<scalar>("vof2IsoTol", 1e-8);
const scalar surfCellTol = isoAdvectorControl.lookupOrDefault<scalar>("surfCellTol", 1e-8);
const scalar clipAlphaTol = isoAdvectorControl.lookupOrDefault<scalar>("clipAlphaTol", 1e-12);
const bool boundAlpha = isoAdvectorControl.lookupOrDefault<bool>("boundAlpha", true);
const bool snapAlpha = isoAdvectorControl.lookupOrDefault<bool>("snapAlpha", true);
const bool writeToLog = isoAdvectorControl.lookupOrDefault<bool>("writeToLog", true);