// Advancing alpah1 field
if (interfaceMethod == "MULES")
{
    #include "alphaEqnMULES.H"
}
else if (interfaceMethod == "isoAdvector")
{
    // Temporarily making U relative to mesh motion
    if (mesh.moving())
    {
        U -= fvc::reconstruct(mesh.phi());
    }

    // Updating alpha1
    advector.advect();

    // Making U absolute again after advection step
    if (mesh.moving())
    {
        U += fvc::reconstruct(mesh.phi());
    }

    const dimensionedScalar& rho1f(rho1);
    const dimensionedScalar& rho2f(rho2);
    rhoPhi = advector.getRhoPhi(rho1f, rho2f);

    alpha2 = 1.0 - alpha1;
    mixture.correct();

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  1 - Max(" << alpha1.name() << ") = " << 1 - max(alpha1).value()
        << endl;
}

//U not solved but prescribed and optionally periodic
if (prescribedU)
{

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
        phi = uFactor*phi.prevIter();
        U = uFactor*U.prevIter();
    }
}
