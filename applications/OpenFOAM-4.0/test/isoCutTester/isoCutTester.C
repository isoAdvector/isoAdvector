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
    isoCutFaceTester

Description
    Testing functions of the isoCutFace class.

Author
    Johan Roenby, DHI, all rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "isoCutFace.H"
#include "isoCutCell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void writePly
(
    const DynamicList< List<point> >& faces,
    const word fileName,
    const word fileDir
)
{
    //Writing faces to ply file for inspection in paraview
    
    autoPtr<OFstream> plyFilePtr;
    plyFilePtr.reset(new OFstream(fileDir + "/" + fileName + ".ply"));
    plyFilePtr() << "ply" << endl;
    plyFilePtr() << "format ascii 1.0" << endl;
    plyFilePtr() << "comment " << fileName << endl;
    label nPoints(0);
    forAll(faces,fi)
    {
        nPoints += faces[fi].size();
    }

    plyFilePtr() << "element vertex " << nPoints << endl;
    plyFilePtr() << "property float32 x" << endl;
    plyFilePtr() << "property float32 y" << endl;
    plyFilePtr() << "property float32 z" << endl;
    plyFilePtr() << "element face " << faces.size() << endl;
    plyFilePtr() << "property list uint8 int32 vertex_index" << endl;
    plyFilePtr() << "end_header" << endl;
    forAll(faces,fi)
    {
        List<point> pf = faces[fi];
        forAll(pf,pi)
        {
            point p = pf[pi];
            plyFilePtr() << p[0] << " " << p[1] << " " << p[2] << endl;
        }
    }
    label np = 0;
    forAll(faces,fi)
    {
        nPoints = faces[fi].size();
        plyFilePtr() << nPoints;
        for (label pi = np; pi < np + nPoints; pi++ )
        {
            plyFilePtr() << " " << pi;
        }
        plyFilePtr() << "" << endl;
        np += nPoints;
    }    
}

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

	const scalarField& x = mesh.points().component(0);
    const scalarField& y = mesh.points().component(1);
    const scalarField& z = mesh.points().component(2);
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

    //Testing isoCutFace
    isoCutFace icf(mesh,f);
    isoCutFace icf2(mesh,f);

    //Subfaces calculated using the fvMesh and scalarField of isoCutFace
    DynamicList< List<point> > subFaces(mesh.nFaces());
    //Subfaces calculated by specifying face points and values in calcSubFace
    DynamicList< List<point> > subFaces2(mesh.nFaces());
    const pointField& points = mesh.points();

    forAll(mesh.faces(),fi)
    {
        label faceStatus = icf.calcSubFace(fi,f0);
        if (faceStatus == 0)
        {
            List<point> sfpi = icf.subFacePoints();
            subFaces.append(sfpi);
            Info << "subfacepoints for face " << fi << ": " << sfpi << endl;
            label nPoints = sfpi.size();
            forAll(sfpi,pi)
            {
                if (mag(sfpi[pi] - sfpi[(pi + 1) % nPoints]) < 1e-10)
                {
                    Info << "Warning: Possible dublicate points for subface "
                        << "of face " << fi << " with owner " 
                        << mesh.owner()[fi] << ". Diff = " 
                        << sfpi[(pi + 1) % nPoints] - sfpi[pi] << endl;
                }
                if (mag(mag(sfpi[pi]-centre) - radius) > .1)
                {
                    Info << "Warning: A point is far from isoface"
                        << " for face " << fi << " with owner " 
                        << mesh.owner()[fi] << ". |r-r0| = " 
                        << mag(mag(sfpi[pi]-centre) - radius) << endl;                    
                }
            }
        }

        const labelList& pLabels = mesh.faces()[fi];
        pointField fPts(pLabels.size());
        scalarField fVals(pLabels.size());
        forAll(pLabels, pi)
        {
            fPts[pi] = points[pLabels[pi]];
            fVals[pi] = f[pLabels[pi]];
        }
        faceStatus = icf.calcSubFace(fPts,fVals,f0);
        if (faceStatus == 0)
        {
            List<point> sfpi = icf.subFacePoints();
            subFaces2.append(sfpi);
            Info << "subfacepoints for face " << fi << ": " << sfpi << endl;
            label nPoints = sfpi.size();
            forAll(sfpi,pi)
            {
                if (mag(sfpi[pi] - sfpi[(pi + 1) % nPoints]) < 1e-10)
                {
                    Info << "Warning: Possible dublicate points for subface "
                        << "of face " << fi << " with owner " 
                        << mesh.owner()[fi] << ". Diff = " 
                        << sfpi[(pi + 1) % nPoints] - sfpi[pi] << endl;
                }
                if (mag(mag(sfpi[pi]-centre) - radius) > .1)
                {
                    Info << "Warning: A point is far from isoface"
                        << " for face " << fi << " with owner " 
                        << mesh.owner()[fi] << ". |r-r0| = " 
                        << mag(mag(sfpi[pi]-centre) - radius) << endl;                    
                }
            }
        }
    }
    
    writePly(subFaces,"subFaces", ".");   
    writePly(subFaces2,"subFaces2",".");    

    //Testing isoCutCell
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
    writePly(isoFaces,"isoFaces",".");    
    
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
