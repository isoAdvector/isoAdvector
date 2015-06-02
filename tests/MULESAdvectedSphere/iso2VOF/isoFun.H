    //Define function on mesh points and isovalue
    const scalarField x = mesh.points().component(0);
    const scalarField y = mesh.points().component(1);
    const scalarField z = mesh.points().component(2);
    const scalarField f = pow(0.2545234,2) - pow(x-.5,2) - pow(y+.6,2) - pow(z-.3,2);
    const scalar f0 = 0;

