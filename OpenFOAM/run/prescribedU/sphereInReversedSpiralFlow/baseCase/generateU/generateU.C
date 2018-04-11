/*---------------------------------------------------------------------------*\
              Original work | Copyright (C) 2016-2017 DHI
              Modified work | Copyright (C) 2016-2017 OpenCFD Ltd.
              Modified work | Copyright (C) 2017-2018 Johan Roenby
-------------------------------------------------------------------------------

License
    This file is part of IsoAdvector, which is an unofficial extension to
    OpenFOAM.

    IsoAdvector is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    IsoAdvector is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with IsoAdvector. If not, see <http://www.gnu.org/licenses/>.

Application
    generateU

Description
    Generates velocity field for the classical test case with a sphere
    deformed into a spiralling sheet and back again.

Author
    Johan Roenby, STROMNING, all rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

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

    const scalar pi = Foam::constant::mathematical::pi;

    const scalarField x(mesh.C().component(vector::X));
    const scalarField y(mesh.C().component(vector::Y));
    const scalarField z(mesh.C().component(vector::Z));
    const scalarField r(sqrt(sqr(x - 0.5) + sqr(y - 0.5)));

    vectorField& Uc = U.primitiveFieldRef();

    Uc.replace(vector::X, sin(2*pi*y)*sqr(sin(pi*x)));
    Uc.replace(vector::Y, -sin(2*pi*x)*sqr(sin(pi*y)));
    Uc.replace(vector::Z, sqr(1.0 - 2.0*r));

    U.correctBoundaryConditions();

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

    // Calculating phi
    const vectorField Cf(mesh.Cf().primitiveField());
    const scalarField Xf(Cf.component(vector::X));
    const scalarField Yf(Cf.component(vector::Y));
    const scalarField Zf(Cf.component(vector::Z));
    const scalarField Rf(sqrt(sqr(Xf - 0.5) + sqr(Yf - 0.5)));

    vectorField Uf(Xf.size());
    Uf.replace(0, sin(2*pi*Yf)*sqr(sin(pi*Xf)));
    Uf.replace(1, -sin(2*pi*Xf)*sqr(sin(pi*Yf)));
    Uf.replace(2, sqr(1.0 - 2.0*Rf));

    scalarField& phic = phi.primitiveFieldRef();
    const vectorField& Sfc = mesh.Sf().primitiveField();
    phic = Uf & Sfc;

    surfaceScalarField::Boundary& phibf = phi.boundaryFieldRef();
    const surfaceVectorField::Boundary& Sfbf =
        mesh.Sf().boundaryField();
    const surfaceVectorField::Boundary& Cfbf =
        mesh.Cf().boundaryField();

    forAll(phibf, patchi)
    {
        scalarField& phif = phibf[patchi];
        const vectorField& Sff = Sfbf[patchi];
        const vectorField& Cff = Cfbf[patchi];
        const scalarField xf(Cff.component(vector::X));
        const scalarField yf(Cff.component(vector::Y));
        const scalarField zf(Cff.component(vector::Z));
        const scalarField rf(sqrt(sqr(xf - 0.5) + sqr(yf - 0.5)));

        vectorField Uf(xf.size());
        Uf.replace(0, sin(2*pi*yf)*sqr(sin(pi*xf)));
        Uf.replace(1, -sin(2*pi*xf)*sqr(sin(pi*yf)));
        Uf.replace(2, sqr(1.0 - 2.0*rf));

        phif = Uf & Sff;
    }

    ISstream::defaultPrecision(18);

    U.write();
    phi.write();

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
