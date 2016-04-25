/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of the IsoAdvector source code library, which is an 
	unofficial extension to OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    isoCutCellTester

Description
    Testing functions of the isoCutCell class.

Author
    Johan Roenby, DHI, all rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "isoCutCell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

	Info<< "Reading isoSurfDict\n" << endl;

	IOdictionary isoSurfDict
	(
		IOobject
		(
			"isoSurfDict",
			runTime.system(),
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	);

	const word surfType = isoSurfDict.lookup("type");
	const vector centre = isoSurfDict.lookup("centre");
	const vector direction = isoSurfDict.lookupOrDefault<vector>("direction",vector::zero);
	const scalar radius = isoSurfDict.lookupOrDefault<scalar>("radius",0.0);

	const scalarField x = mesh.points().component(0);
    const scalarField y = mesh.points().component(1);
    const scalarField z = mesh.points().component(2);
    scalar f0 = 0.0;
	scalarField f(x.size());

	if ( surfType == "plane" )
	{
		f = -(mesh.points() - centre) & (direction/mag(direction));
		f0 = 0.0;
	}
	else if ( surfType == "sphere" )
	{
		f = Foam::exp(-pow((x-centre[0]),2) - pow((y-centre[1]),2) - pow((z-centre[2]),2));
		f0 = Foam::exp(-pow(radius,2));
	}
	else if ( surfType == "cylinder" )
	{
		f = -sqrt(pow(mag(mesh.points()-centre),2) - pow(mag((mesh.points()-centre) & direction),2));
		f0 = -radius;
	}
	else
	{
		Info << "Invalid surface type specified" << endl;
		Info << "Aborting..." << endl;
	}

    //Define function on mesh points and isovalue
	    
    isoCutCell icc(mesh,f);
        
    DynamicList< List<point> > isoFaces(mesh.nCells());
    forAll(mesh.cells(),ci)
    {
        label cellStatus = icc.calcSubCell(ci,f0);
        Info << "Cell status for cell " << ci << " is " << cellStatus << endl;
        if (cellStatus == 0)
        {
            List<point> ifpi = icc.isoFacePoints();
            Info << "Cell " << ci << " is cut with isoFace centre: "
                << icc.isoFaceCentre() << " and isoFace area: "
                << icc.isoFaceArea() << " and isoFace points: " << ifpi
                << endl;
            isoFaces.append(ifpi);
            label nPoints = ifpi.size();
            forAll(ifpi,pi)
            {
                if (mag(ifpi[pi] - ifpi[(pi + 1) % nPoints]) < 1e-10)
                {
                    Info << "Warning: Possible dublicate points for isoface "
                        << " of cell " << ci << ": " << ifpi[pi] << " and " 
                        << ifpi[(pi + 1) % nPoints] << endl;
                }
            }
        }
    }
    
    //Writing isofaces to ply file for inspection in paraview
    word fileName = "isoFaces";
    word fileDir = ".";
    
    autoPtr<OFstream> plyFilePtr;
    plyFilePtr.reset(new OFstream(fileDir + "/" + fileName + ".ply"));
    plyFilePtr() << "ply" << endl;
    plyFilePtr() << "format ascii 1.0" << endl;
    plyFilePtr() << "comment " << fileName << endl;
    label nPoints(0);
    forAll(isoFaces,fi)
    {
        nPoints += isoFaces[fi].size();
    }

    plyFilePtr() << "element vertex " << nPoints << endl;
    plyFilePtr() << "property float32 x" << endl;
    plyFilePtr() << "property float32 y" << endl;
    plyFilePtr() << "property float32 z" << endl;
    plyFilePtr() << "element face " << isoFaces.size() << endl;
    plyFilePtr() << "property list uint8 int32 vertex_index" << endl;
    plyFilePtr() << "end_header" << endl;
    forAll(isoFaces,fi)
    {
        List<point> pf = isoFaces[fi];
        forAll(pf,pi)
        {
            point p = pf[pi];
            plyFilePtr() << p[0] << " " << p[1] << " " << p[2] << endl;
        }
    }
    label np = 0;
    forAll(isoFaces,fi)
    {
        nPoints = isoFaces[fi].size();
        plyFilePtr() << nPoints;
        for (label pi = np; pi < np + nPoints; pi++ )
        {
            plyFilePtr() << " " << pi;
        }
        plyFilePtr() << "" << endl;
        np += nPoints;
    }
    
	
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
