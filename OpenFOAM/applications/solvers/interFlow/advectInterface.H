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
}