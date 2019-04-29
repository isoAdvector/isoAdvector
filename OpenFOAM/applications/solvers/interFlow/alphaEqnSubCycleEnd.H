
if (interfaceMethod == "isoAdvector" && !pimple.firstPimpleIter())
{
    // Resetting U and phi to value at latest iteration.
    U = 2.0*U - U.prevIter();
    phi = 2.0*phi - phi.prevIter();
}

