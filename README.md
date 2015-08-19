# IsoAdvection Project  

This is the repository of the IsoAdvection project. 

Developed by Johan Rønby <jro@dhigroup.com>  

Contributors:  

Tomislav Maric <tomislav@sourceflux.de> 

# Project structure 

`src/` 

* Contains a single library : `libIsoAdvection`. 
* Two classes implement the advection scheme: `isoAdvection` and `isoCutter`. 

`applications/` 

- `applications/utilities/postProcessing/isoCut` 
    - Sets the initial volume fraction field. 
- `applications/solvers/isoAdvect` 
    - Solves the volume fraction advection equation. 

## IsoAdvection 

## IsoCutter

- This OpenFOAM code calculates an isosurface from a function f and a function value f0.
- The function is defined in all vertex points of an fvMesh.
- The routine travels throug all cells and looks at each cell face edge.
- If the f is above f0 at one vertex and below f0 at the other vertex of an edge, the edge is cut at a position determined by linear interpolation.
- The routine calculates the polygonal "isoFace" inside the cell formed by cutting its edges in this way.
- It also calculates the volume and cell centre of the "submerged" subcell defined as the part of the original cell where f >= f0.

So far the routine assumes that a cut face is cut at two edges. For n-gons with n > 3 there will in general be special situations where this is not the case. The function should then triangulate such faces interpolating the function to the face centre and do the face cutting on the triangular faces.

The routine was deliberately build to travel through all cells instead of all faces or edges to obtain its results. The downside of this strategy is that the number of times an edge is visited is equal to twice the number of faces to which it belongs.  The advantage is that parallelisation is trivial. It should be noted that the operation performed on each edge is extremely simple.

# Usage  

1. `isoCut`
2. `isoAdvect`
3. ... 

## Compilation 

1. Source the `etc/bashrc` script to set the `ISOADVECTION` project environment variable. 

~~~
    source etc/bashrc
~~~

2. Execute the `Allwmake` script. 

## Preprocessing 

## Running 

## Postprocessing 

# Ongoing work 

Add feature and improvement ideas/requests here. Some of them will be taken over by the development team and converted to issues. If you see an open feature request in this list that you want to take on, open an issue and leave a comment in the issue tracker.  

## Bring the code to the level before the HD failure:

- Setup notchedDisk/Zalesak test case (Remember stream function approach to calculating fluxes from analytical velocity field)
- Setup scalarTransportFoam cases with HRIC, CICSAM, upwind, vanLeer01 (get setups from navalFoam cases)

## To bring the code further:

- Setup unstructured version of advectedSphere
- Reimplement method to find isovalue that gives VOF value
- Make switch that allows one to turn on/off debug messages from isoAdvector

- Implement isoAdvector in interFoam instead of MULES: isoInterFoam
- Setup filling tank case
- Setup standing wave test case 

- Develop and implement consistent boundary conditions for isosurface reconstruction step.

### `isoCutter` enhancements 

- Profile isocutter with valgrind
- Think about how to replace all the DynamicList if necessary
- Add function to find isovalue that returns a given vof value

- Arrange data and functions into classes (done)
- Remove terminal write out of face points (maybe add point fields to write out for paraview) (done)
- Make behaviour correct in special cases where vertices are exactly on the isoSurface (done)
- Make a .ply polygon face writer for viewing iso surface in paraview (done)
- Make function write a volScalarField for visualisation with paraview (done)
- Add function that calculates submerged volume for each cell (done)
	- Add function calculating face centres and areas from new face points (done)
	- Add function that extracts face centres and areas from fully submerged faces (done)
	- Add function that calculates volume from this taking care of area vector direction (done)
- Add function to find fully submerged faces of each cut cell (done)
- Add function to calculate submerged part of cut faces for each cell (done)
- Add indicator of orientation is correct (done - and undone)
- make an isoCutCell(cellI) instead of the current isoCutCells (done)

