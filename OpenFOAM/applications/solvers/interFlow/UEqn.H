// Looking up which alpha1 solver to use
const dictionary& UDict(mesh.solverDict(U.name()));

bool prescribedU
(
    UDict.lookupOrDefault<bool>("prescribedU", false)
);

//U not solved but prescribed and optionally periodic
if (prescribedU)
{
    const scalar period(UDict.lookupOrDefault<scalar>("period", 0));
    scalar reverseTime(UDict.lookupOrDefault<scalar>("reverseTime", 0));

    //Setting velocity field and face fluxes for next time step
    const scalar t(runTime.time().value());
    const scalar dt(runTime.deltaT().value());
    if (reverseTime > 0.0 && t >= reverseTime)
    {
        Info<< "Reversing flow" << endl;
        phi *= -1.0;
        U *= -1.0;
        reverseTime = -1.0;
    }
    if (period > 0.0)
    {
        const scalar uFactor
        (
            .5*(Foam::cos(2.0*M_PI*t/period)
                + Foam::cos(2.0*M_PI*(t + dt)/period))
        );
        phi = uFactor*phi0;
        U = uFactor*U0;
    }
}
else
{
    #include "UEqnSolved.H"
}