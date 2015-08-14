    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"

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
	const scalar x0(0.5), z0(0.5), rad(0.25);
    const scalarField x = mesh.points().component(0);
    const scalarField y = mesh.points().component(1);
    const scalarField z = mesh.points().component(2);
	const scalar Ux = U[0].component(0); //Assumes spatially constant U
	const scalar Uz = U[0].component(2); //Assumes spatially constant U
	scalar xc = x0 + Ux*runTime.time().value();
	scalar zc = z0 + Uz*runTime.time().value();
    scalarField f = -sqrt(pow(x-xc,2) + pow(z-zc,2));;
	
	//Calculating alpha1 volScalarField from f = f0 isosurface
    Foam::isoCutter cutter(mesh,f);
    cutter.subCellFractions(-rad,alpha1);
	alpha1.correctBoundaryConditions();
    alpha1.write(); //Writing volScalarField alpha1

	Info << "sum(alpha*V) = " << sum(mesh.V()*alpha1).value() 
	 << ", max(alpha1)-1 = " << max(alpha1).value()-1.0
	 << "\t min(alpha1) = " << min(alpha1).value() << endl;