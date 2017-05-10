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
        phi = uFactor*phi0;
        U = uFactor*U0;
    }
}

if (interfaceMethod == "isoAdvector")
{
    #include "isoAdvectInterface.H"
}
else if (interfaceMethod == "MULES")
{
    #include "alphaEqnSubCycle.H"
}
else if (interfaceMethod == "fvSchemes")
{
    solve
    (
        fvm::ddt(alpha1)
      + fvm::div(phi, alpha1)
     ==
        fvOptions(alpha1)
    );

    //Clip and snap alpha1 to ensure strict boundedness to machine precision
    alpha1Unbounded = alpha1;
    if ( snapAlphaTol > 0.0 )
    {
        alpha1 = alpha1*
            pos(alpha1-snapAlphaTol)*neg(alpha1-(1.0-snapAlphaTol))
            + pos(alpha1-(1.0-snapAlphaTol));
    }
    if ( clipAlpha )
    {
        alpha1 = min(1.0,max(0.0,alpha1));
    }

    //Calculate rho and rhoPhi
    rho == alpha1*rho1 + (scalar(1) - alpha1)*rho2;
    rhoPhi = advector.getRhoPhi(rho1, rho2);

}
