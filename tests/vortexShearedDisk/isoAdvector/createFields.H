	Info<< "Reading field alpha1\n" << endl;
    volScalarField alpha1
    (
        IOobject
        (
            "alpha1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    //Define function on mesh points and isovalue
    const scalarField x = mesh.points().component(0);
    const scalarField y = mesh.points().component(1);
    const scalarField z = mesh.points().component(2);
    const scalarField f = pow(0.15,2) - pow(x-.5,2) - pow(z-.75,2);
    const scalar f0 = 0;
	
	//Calculating alpha1 volScalarField from f = f0 isosurface
    Foam::isoCutter cutter(mesh,f);
    cutter.subCellFractions(f0,alpha1);
	alpha1.correctBoundaryConditions();
    alpha1.write(); //Writing volScalarField alpha1

	Info << "sum(alpha*V) = " << sum(mesh.V()*alpha1).value() 
	 << ", max(alpha1)-1 = " << max(alpha1).value()-1.0
	 << "\t min(alpha1) = " << min(alpha1).value() << endl;
	 
    Info<< "Reading field U\n" << endl;
		
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

	scalar PI = Foam::constant::mathematical::pi;
	scalarField X = mesh.C().component(0);
	scalarField Z = mesh.C().component(2);
	scalarField u = sin(2*PI*Z)*pow(sin(PI*X),2);
	scalarField w = -sin(2*PI*X)*pow(sin(PI*Z),2);
	forAll(U,ci)
	{
		U[ci] = u[ci]*vector(1,0,0) + w[ci]*vector(0,0,1);
	}
	
	Info<< "Reading/calculating face flux field phi\n" << endl;

	surfaceScalarField phi
	(
		IOobject
		(
			"phi",
			runTime.timeName(),
			mesh,
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
		linearInterpolate(U) & mesh.Sf()
	);

	scalarField psi = -pow(sin(PI*x)*sin(PI*z),2)/PI;
	forAll(phi,fi)
	{
		phi[fi] = 0.0;
		labelList faceLabels = mesh.faces()[fi];
		label nPoints = faceLabels.size();
		forAll(faceLabels,pi)
		{
			label pl1 = faceLabels[pi];
			label pl2 = faceLabels[(pi+1) % nPoints];
			phi[fi] += 0.5*(psi[pl1]+psi[pl2])*(y[pl2]-y[pl1]);
		}
	}
	
	U.write();
	phi.write();