// Instantiating isoAdvection object for advancing alpha1
isoAdvection advector(alpha1, phi, U);

// Checking if U is to be treated as prescribed
const dictionary& UDict(mesh.solverDict(U.name()));

bool prescribedU
(
    UDict.lookupOrDefault<bool>("prescribedU", false)
);

// Checking if U is to be treated as periodic or reversing
const scalar period(UDict.lookupOrDefault<scalar>("period", 0));
scalar reverseTime(UDict.lookupOrDefault<scalar>("reverseTime", 0));

//Creating copies of initial U and phi for periodic flows
//Rewrite to  avoid these extra fields when predcribedU = false
if (prescribedU)
{
//const volVectorField U0(U);
//const surfaceScalarField phi0(phi);
    U.storePrevIter();
    phi.storePrevIter();
}
