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

advectionTime += (runTime.elapsedCpuTime() - advectionStartTime);

Info<< "Advection time is " << 100*advectionTime/runTime.elapsedCpuTime() 
    << "% of elapsedCpuTime." << endl;

const scalar V = gSum(mesh.V().field()*alpha1.internalField());
Info<< "Total alpha1 volume: " << V << " m3. Deviation from initial: " 
    << 100*(V-V0)/V0 << "%." << nl << endl;