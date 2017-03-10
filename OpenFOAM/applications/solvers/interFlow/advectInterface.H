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
