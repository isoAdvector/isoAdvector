# IsoAdvection Project  

This is the repository of the IsoAdvection project. 

Developed by Johan Roenby <jro@dhigroup.com>  

Contributors:

* Hrvoje Jasak <hrvoje.jasak@fsb.hr> (Code clean up, consistent treatment of boundary faces including processor boundaries)
* Vuko Vukcevic <vuko.vukcevic@fsb.hr> (Code review and profiling)
* Tomislav Maric <tomislav@sourceflux.de> (Source file rearrangement)

# Project structure 

`src/` 

* Contains a single library : `libIsoAdvection`. 
* Two classes implement the advection scheme: `isoAdvection` and `isoCutter`. 

`applications/` 

- `applications/solvers/isoAdvect` 
    - Solves the volume fraction advection equation in either steady flow or periodic flow with the option of changing the flow direction at a specified time.
- `applications/utilities/preProcessing/isoSurf` 
    - Sets the initial volume fraction field for either a sphere, a cylinder or a plane. 
- `applications/utilities/postProcessing/uniFlowErrors`
    - For cases with spheres and discs in steady uniform flow calculates errors relative to exact VOF solution

## IsoAdvection 

- Calculates the total volume of water crossing each face in the mesh during the time interval from time t to time t + dt.

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

See the tests directory for examples of usage.
For instance go to tests/discInUniFlow/hex and execute the generateCase script to make a simple case.
Go to the newly created case folder and execute the Allrun script.
After a few seconds it should finish and inspection of the alpha1 field in Paraview will show a circular water volume going from lower left to upper right corner of the domain.

## Compilation 

Execute the `Allwmake` script. 

To write (A LOT OF!) debug data to log file write "./Allwmake debug" when compiling.

## Preprocessing 

## Running 

## Postprocessing 

# Ongoing work 

Add feature and improvement ideas/requests here. Some of them will be taken over by the development team and converted to issues. If you see an open feature request in this list that you want to take on, open an issue and leave a comment in the issue tracker.  

## To bring the code further:

- Reimplement method to find isovalue that gives VOF value
- Setup notchedDisk/Zalesak test case (Remember stream function approach to calculating fluxes from analytical velocity field)
- Setup filling tank case
- Setup standing wave test case 
