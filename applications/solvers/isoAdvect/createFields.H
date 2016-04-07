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

	//Volume of water transported across faces in time dt
	surfaceScalarField dVf
	(
		IOobject
		(
			"dVf",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("vol", dimVol, 0)
    );	

	//Initial Volume of Fluid (VOF)
    scalar V0 = sum(mesh.V()*alpha1).value();   

    Info << "sum(alpha*V) = " << V0
     << ", max(alpha1)-1 = " << max(alpha1).value()-1.0
     << "\t min(alpha1) = " << min(alpha1).value() << endl;

    const dictionary& isoAdvectorDict = mesh.solutionDict().subDict("isoAdvector");
    const scalar clipAlphaTol = isoAdvectorDict.lookupOrDefault<scalar>("clipAlphaTol", 1e-12);
    const bool snapAlpha = isoAdvectorDict.lookupOrDefault<bool>("snapAlpha", true);
	const scalar period = isoAdvectorDict.lookupOrDefault<scalar>("period", 0.0);
	scalar reverseTime = isoAdvectorDict.lookupOrDefault<scalar>("reverseTime", 0.0);

	//Copies of initial U and phi for use when flow is periodic
	surfaceScalarField phi0 = phi;
	volVectorField U0 = U;
	scalar PI = Foam::constant::mathematical::pi;
