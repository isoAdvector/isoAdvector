Welcome to the IsoAdvector project

# What is IsoAdvector?

IsoAdvector is a geometric Volume-of-Fluid method for advection of a sharp interface between two incompressible fluids.
It works on both structured and unstructured meshes with no requirements on cell shapes.
IsoAdvector is meant as a replacement for the MULES scheme implemented in the interFoam family of solvers.

# Contributors:

* Johan Roenby <jro@dhigroup.com> (Inventor and main developer)
* Hrvoje Jasak <hrvoje.jasak@fsb.hr> (Consistent treatment of boundary faces including processor boundaries, code clean up)
* Vuko Vukcevic <vuko.vukcevic@fsb.hr> (Code review and profiling)
* Tomislav Maric <tomislav@sourceflux.de> (Source file rearrangement)

# Requirement:

IsoAdvector is implemented with OpenFOAM-2.2.0. 
It will also work with other OpenFOAM/foam-extend versions with minor adjustments by the user.

# Installation:

1.  In a linux terminal download the package with git by typing:

        git clone https://bitbucket.org/dhifoam/isoadvector.git

2.  Execute the Allwmake script by typing (will finish in a ~1 min):

        ./Allwmake

    (To write A LOT OF debug data do ./Allwmake debug)
	
3.  Test the installation with a simple test case by typing (will finish in seconds):
    
	    cd tests/discInUniFlow
		./makeAndRunTestCase
	
    Open Paraview and color the mesh by alpha1.

# Code structure:

`src/` 

* Contains a single library : `libIsoAdvection`. 
* Two classes implement the IsoAdvector advection scheme: `isoAdvection` and `isoCutter`. 
* For comparison the library includes HRIC, CICSAM and other algebraic VOF scheme in `finiteVolume` directory. 

`applications/` 

- `applications/solvers/isoAdvect` 
    - Solves the volume fraction advection equation in either steady flow or periodic predefined flow with the option of changing the flow direction at a specified time.
- `applications/solvers/mulesFoam`
    - This is like isoAdvect, but using MULES instead of isoAdvector. Based on interFoam.
- `applications/solvers/passiveAdvectionFoam` 
    - This is essentially scalarTransportFoam but using alpha1 instead of T and without diffusion term. Used to run test cases with predefined velocity field using the CICSAM and HRIC schemes.
- `applications/solvers/interFlow` 
    - This is essentially the original interFoam solver with MULES replaced by isoAdvector for interface advection. No changes to PIMPLE loop.
- `applications/utilities/preProcessing/isoSurf` 
    - Sets the initial volume fraction field for either a sphere, a cylinder or a plane. See isoSurfDict for usage.
- `applications/utilities/postProcessing/uniFlowErrors`
    - For cases with spheres and discs in steady uniform flow calculates errors relative to exact VOF solution.

`applications/`

* Contains a series of test cases
	
#Classes
	
## IsoAdvection 

- Calculates the total volume of water, dVf, crossing each face in the mesh during the time interval from time t to time t + dt.

## IsoCutter

- Calculates an isosurface from a function f and a function value f0.
- The function is defined in all vertex points of an fvMesh.
- The routine travels throug all cells and looks at each cell face edge.
- If the f is above f0 at one vertex and below f0 at the other vertex of an edge, the edge is cut at a position determined by linear interpolation.
- The routine calculates the polygonal "isoFace" inside the cell formed by cutting its edges in this way.
- It also calculates the volume and cell centre of the "submerged" subcell defined as the part of the original cell where f >= f0.

So far the routine assumes that a cut face is cut at two edges. For n-gons with n > 3 there will in general be special situations where this is not the case. The function should then triangulate such faces interpolating the function to the face centre and do the face cutting on the triangular faces.

The routine was deliberately build to travel through all cells instead of all faces or edges to obtain its results. The downside of this strategy is that the number of times an edge is visited is equal to twice the number of faces to which it belongs.  The advantage is that parallelisation is trivial. It should be noted that the operation performed on each edge is extremely simple.

# Ongoing work 

Add feature and improvement ideas/requests here. Some of them will be taken over by the development team and converted to issues. If you see an open feature request in this list that you want to take on, open an issue and leave a comment in the issue tracker.  

## To bring the code further:

- Reimplement method to find isovalue that gives VOF value
- Setup notchedDisk/Zalesak test case (Remember stream function approach to calculating fluxes from analytical velocity field)
- Setup filling tank case
- Setup standing wave test case 
