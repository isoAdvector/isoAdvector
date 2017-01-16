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

//Advance alpha1 from time t to t+dt
scalar advectionStartTime = runTime.elapsedCpuTime();

advector.advect();

//Print surface cell set and bounded cell set if required by user
if (printSurfCells)
{
    advector.getSurfaceCells(surfCells);
}

if (printBoundCells)
{
    advector.getBoundedCells(boundCells);
}

//Clip and snap alpha1 to ensure strict boundedness to machine precision
alpha1Unbounded = alpha1;
if ( clipAlphaTol > 0.0 )
{
    alpha1 = alpha1*
        pos(alpha1-clipAlphaTol)*neg(alpha1-(1.0-clipAlphaTol)) 
        + pos(alpha1-(1.0-clipAlphaTol));
}
if ( snapAlpha )
{
    alpha1 = min(1.0,max(0.0,alpha1));
}

//Calculate rho and rhoPhi
rho == alpha1*rho1 + (scalar(1) - alpha1)*rho2;
rhoPhi = advector.getRhoPhi(rho1, rho2);

advectionTime += (runTime.elapsedCpuTime() - advectionStartTime);

Info<< "Advection time is " << 100*advectionTime/runTime.elapsedCpuTime() 
    << "% of elapsedCpuTime." << endl;

const scalar V = gSum(mesh.V().field()*alpha1.internalField());
Info<< "Total alpha1 volume: " << V << " m3. Deviation from initial: " 
    << 100*(V-V0)/V0 << "%." << nl << endl;