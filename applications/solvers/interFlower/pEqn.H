{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("Dp", fvc::interpolate(rAU));

    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn.H();

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (fvc::interpolate(HbyA) & mesh.Sf())
      + fvc::ddtPhiCorr(rAU, rho, U, phi)
    );

    adjustPhi(phiHbyA, U, p_rgh);
    phi = phiHbyA;
/*
	Info << "Creating xSnSdotnfdl" << endl;
	surfaceVectorField xSnSdotnfdl
	(
		IOobject
		(
			"xSnSdotnfdl",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
	    dimensionedVector("area", dimArea, vector::zero)
	);
	
	Info << "Calculating xSnSdotnfdl" << endl;
	isoAdvection advector2(alpha1,phi,U,isoAdvectorDict);
	advector2.getxSnSdotnF(xSnSdotnfdl);

	
	Info << "Creating Ci" << endl;
	volVectorField Ci
	(
		IOobject
		(
			"Ci",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
	    dimensionedVector("len", dimLength, vector::zero)
	);
	Info << "Creating Si" << endl;
	volVectorField Si
	(
		IOobject
		(
			"Si",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
	    dimensionedVector("area", dimArea, vector::zero)
	);
	
	advector2.getIsoCentreAndNormal(Ci,Si);
	
	Info << "max(mag(xSnSdotnfdl)) = " << max(mag(xSnSdotnfdl)) << endl;
	Info << "max(mag(Ci)) = " << max(mag(Ci)) << endl;
	Info << "max(mag(Si)) = " << max(mag(Si)) << endl;
	
	Info << "Using	 xSnSdotnfdl" << endl;
*/
//	volVectorField Ug('Ug',((Ci & g)*((rho1-rho2)/mesh.V())*Si));
//	volVectorField Ug = (Ci & g)*Si;
//	surfaceScalarField phig2 = ((Ci & g)*(rho1-rho2))*Si;
    surfaceScalarField phig
    (
        (
            fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
//            fvc::interpolate(rAU*interface.sigmaK())*fvc::snGrad(alpha1)*mesh.magSf()
//          - (fvc::interpolate(gh*rAU*fvc::grad(rho)) & mesh.Sf())
          - ghf*fvc::snGrad(rho)
//          - fvc::snGrad(gh*rho)
        )*rAUf*mesh.magSf()
//        - (rho2-rho1)*(g & xSnSdotnfdl)*rAUf
//		- (rho1-rho2)*(fvc::interpolate(Ug) & mesh.Sf())*rAUf
//        )*mesh.magSf()
//        )
    );

/*	
	surfaceScalarField v1
	(
		IOobject
		(
			"v1",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
	    - ghf*fvc::snGrad(rho)*mesh.magSf()
	);
	v1.write();

	surfaceScalarField v2
	(
		IOobject
		(
			"v3",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
	    (rho2-rho1)*(g & xSnSdotnfdl)
	);
	v2.write();
*/
	
    phiHbyA += phig;

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p_rghEqn.flux();

            U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    #include "continuityErrs.H"

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }
}
