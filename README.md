Welcome to the IsoAdvector project

# What is IsoAdvector?

IsoAdvector is a geometric Volume-of-Fluid method for advection of a sharp 
interface between two incompressible fluids. It works on both structured and 
unstructured meshes with no requirements on cell shapes. IsoAdvector is meant as 
a replacement for the interface compression with the MULES limiter implemented 
in the interFoam family of solvers.

The ideas behind and performance of the IsoAdvector scheme is documented in this
article:

https://arxiv.org/abs/1601.05392

A number of videos can be found in this youtube channel:

https://www.youtube.com/channel/UCt6Idpv4C8TTgz1iUX0prAA


# Requirements:

The isoAdvector library is compatible with several different versions of 
OpenFOAM/foam-extend. The isoAdvector root directory contains a folder named 
after each OpenFOAM/foam-extend version with which the isoAdvcetor code has been 
succesfully compiled. The code will compile with other versions with minor 
modifications. I will do my best to correct bugs in all versions of the code, 
but it is unrealistic that I will have time to implement all changes and new 
developments to older versions of the code. The code is being maintained and 
developed in the newest OpenFOAM version, so as a rule of thumb the best 
performance should be expected with this version, and slight difference in 
behaviour should be expected in older version of the code.

# Installation:

0.  Source a supported OpenFOAM environment: 

        source OpenFOAM-x.x.x/etc/bashrc

1.  In a linux terminal download the package with git by typing:

        git clone https://bitbucket.org/dhifoam/isoadvector.git

2.  Execute the Allwmake script by typing (will finish in a ~1 min):

        cd isoadvector
        ./Allwmake

    Applications will be compiled to your FOAM_USER_APPBIN and libraries will be
    compiled to your FOAM_USER_LIBBIN.
    
3.  Test installation with a simple test case by typing (finishes in secs):
    
	    cp -r [foam version]/run/isoAdvector/discInUniFlow/baseCase ~
        cd ~/baseCase
        ./Allrun
	
    Here [foam version] is the loaded foam version, e.g. OpenFOAM-4.0.
    Open Paraview and color the mesh by the alpha.water/alpha1 field.

If you want to test the code on unstructured meshes, a number of such meshes can
be donwloaded by using the downloadMeshes script in the bin directory. These 
meshes will take up ~2GB and will be placed in a folder called meshes. The 
meshes will be loaded by the scripts in discInUniFlows, vortexSmearedDisc, 
sphereInUniFlow and smearedSphere cases in run/isoAdvector.
    
# Code structure:

`src/` 

* Contains the three classes implement the IsoAdvector advection scheme:
    - `isoCutFace`
    - `isoCutCell` 
    - `isoAdvection` 
  These are compiled into a library named `libIsoAdvection`. 
* For comparison we also include the CICSAM, HRIC and mHRIC algebraic VOF 
  schemes in `finiteVolume` directory. These are compiled into a library called
  libVOFInterpolationSchemes.

`applications/` 

- `applications/solvers/isoAdvector` 
    - Solves the volume fraction advection equation in either steady flow or 
      periodic predefined flow with the option of changing the flow direction at
      a specified time.
- `applications/solvers/mulesFoam`
    - This is like isoAdvector, but using MULES instead of isoAdvector. Based on
      interFoam.
- `applications/solvers/passiveAdvectionFoam` 
    - This is essentially scalarTransportFoam but using alpha1 instead of T and 
      without diffusion term. Used to run test cases with predefined velocity 
      field using the CICSAM and HRIC schemes.
- `applications/solvers/interFlow` 
    - This is essentially the original interFoam solver with MULES replaced by 
      isoAdvector for interface advection. No changes to PIMPLE loop.
- `applications/utilities/preProcessing/isoSurf` 
    - Sets the initial volume fraction field for either a sphere, a cylinder or 
      a plane. See isoSurfDict for usage.
- `applications/utilities/postProcessing/uniFlowErrors`
    - For cases with spheres and discs in steady uniform flow calculates errors 
      relative to exact VOF solution.

`run/`

- `isoAdvector/` 
    - Contains test cases using isoAdvector to move a volume of fluid in a 
      predescribed velocity field.
- `interFlow/` 
    - Contains test cases using interFlow coupling IsoAdvector with the PIMPLE 
      algorithm for the pressure-velocity coupling.
	
#Classes
	
## IsoAdvection 

- Calculates the total volume of water, dVf, crossing each face in the mesh 
  during the time interval from time t to time t + dt.

## IsoCutCell

- Performs cutting of a cell given an isovalue and the alpha values interpolated 
  to the cell vertices.
- Calculates the submerged volume of a cell for a given isovalue.
- Calculates the isovalue given a specified target volume fraction of a cell.
  
## IsoCutFace

- Performs cutting of a face given an isovalue and alpha values interpolated to
  the face vertices.
- Calculates the submerged face area given an isovalue and alpha values 
  interpolated to the face vertices.

So far the routine assumes that a cut face is cut at two edges. For n-gons with 
n > 3 there will in general be special situations where this is not the case. 
The function should then triangulate such faces interpolating the function to 
the face centre and do the face cutting on the triangular faces.

The routine was deliberately build to travel through all cells instead of all 
faces or edges to obtain its results. The downside of this strategy is that the 
number of times an edge is visited is equal to twice the number of faces to 
which it belongs. The advantage is that parallelisation is trivial. It should be 
noted that the operation performed on each edge is extremely simple.


# Contributors:

* Johan Roenby <jro@dhigroup.com> (Inventor and main developer)
* Hrvoje Jasak <hrvoje.jasak@fsb.hr> (Consistent treatment of boundary faces 
  including processor boundaries, code clean up)
* Vuko Vukcevic <vuko.vukcevic@fsb.hr> (Code review, profiling, porting to 
  foam-extend, bug fixing)
* Tomislav Maric <tomislav@sourceflux.de> (Source file rearrangement)
