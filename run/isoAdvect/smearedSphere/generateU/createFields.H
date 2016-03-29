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

	if ( runTime.time().value() <= 0.0 ) //Ugly way to be able to continue a stopped simulation
	{
		//Define function on mesh points and isovalue
		const scalarField x = mesh.points().component(0);
		const scalarField y = mesh.points().component(1);
		const scalarField z = mesh.points().component(2);
		const scalarField f = -sqrt(pow(x-.35,2) + pow(y-.35,2) + pow(z-.35,2));
		const scalar f0 = -0.15;
		
		//Calculating alpha1 volScalarField from f = f0 isosurface
		Foam::isoCutter cutter(mesh,f);
		cutter.subCellFractions(f0,alpha1);
		alpha1.correctBoundaryConditions();
		alpha1.write(); //Writing volScalarField alpha1
	}
	
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
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh
    );

	scalar PI = Foam::constant::mathematical::pi;
	{
		scalarField X = mesh.C().component(0);
		scalarField Y = mesh.C().component(1);
		scalarField Z = mesh.C().component(2);
		scalarField u = 2.0*pow(sin(PI*X),2)*sin(2.0*PI*Y)*sin(2.0*PI*Z);
		scalarField v = -sin(2.0*PI*X)*pow(sin(PI*Y),2)*sin(2.0*PI*Z);
		scalarField w = -sin(2.0*PI*X)*sin(2.0*PI*Y)*pow(sin(PI*Z),2);
		forAll(U,ci)
		{
			U[ci] = u[ci]*vector(1.0,0.0,0.0) + v[ci]*vector(0.0,1.0,0.0) + w[ci]*vector(0.0,0.0,1.0);
		}
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
			IOobject::NO_WRITE
		),
		linearInterpolate(U) & mesh.Sf()
	);

	{
		scalarField Xf = mesh.Cf().component(0);
		scalarField Yf = mesh.Cf().component(1);
		scalarField Zf = mesh.Cf().component(2);
		scalarField uf = 2.0*pow(sin(PI*Xf),2)*sin(2.0*PI*Yf)*sin(2.0*PI*Zf);
		scalarField vf = -sin(2.0*PI*Xf)*pow(sin(PI*Yf),2)*sin(2.0*PI*Zf);
		scalarField wf = -sin(2.0*PI*Xf)*sin(2.0*PI*Yf)*pow(sin(PI*Zf),2);
		forAll(phi,fi)
		{
			vector Uf = uf[fi]*vector(1.0,0.0,0.0) + vf[fi]*vector(0.0,1.0,0.0) + wf[fi]*vector(0.0,0.0,1.0);
			phi[fi] = (Uf & (mesh.Sf()[fi]));
		}
	}

/*
	scalarField psi = (1.0/PI)*pow(sin(PI*x),2.0)*pow(sin(PI*z),2.0);
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
//			phi[fi] += 0.5*(psi[pl1]*y[pl2]-psi[pl2]*y[pl1]);
		}
	}
*/	
	scalarField sumPhi = fvc::surfaceIntegrate(phi);	
	scalar maxMagSumPhi(0.0);
	label maxLabel(0);
	forAll(sumPhi,ci)
	{
		scalar msp = mag(sumPhi[ci]);
		if (msp > maxMagSumPhi )
		{
			maxMagSumPhi = msp;
			maxLabel = ci;
		}
	}
	Info << "maxMagSumPhi/cellVol = " << maxMagSumPhi/mesh.V()[maxLabel] << endl;
	
	U.write();
	phi.write();