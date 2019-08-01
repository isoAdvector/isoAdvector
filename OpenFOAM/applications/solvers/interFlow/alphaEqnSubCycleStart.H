// Looking up which alpha1 solver to use
const dictionary& alpha1Dict(mesh.solverDict(alpha1.name()));

word interfaceMethod
(
    alpha1Dict.lookupOrDefault<word>("interfaceMethod", "isoAdvector")
);

if (interfaceMethod == "isoAdvector" && pimple.nCorrPIMPLE() > 1)
{
    if (pimple.firstPimpleIter())
    {
        // To recalculate the alpha1 update in subsequent iterations, we
        // must store its current value before overwriting with the new
        // value
        alpha1.oldTime();
    }
    else
    {
        // Resetting alpha1 to value before advection in first PIMPLE
        // iteration.
        alpha1 = alpha1.oldTime();
    }
}

