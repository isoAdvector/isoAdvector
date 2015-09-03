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

    Info << "sum(alpha*V) = " << sum(mesh.V()*alpha1).value() 
     << ", max(alpha1)-1 = " << max(alpha1).value()-1.0
     << "\t min(alpha1) = " << min(alpha1).value() << endl;

    const dictionary& isoAdvectorDict = mesh.solutionDict().subDict("isoAdvector");

    const scalar clipAlphaTol = isoAdvectorDict.lookupOrDefault<scalar>("clipAlphaTol", 1e-12);

    const bool snapAlpha = isoAdvectorDict.lookupOrDefault<bool>("snapAlpha", true);