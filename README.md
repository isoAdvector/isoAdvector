# Welcome to the IsoAdvector project

## What is IsoAdvector?

IsoAdvector is a geometric Volume-of-Fluid method for advection of a sharp 
interface between two incompressible fluids. It works on both structured and 
unstructured meshes with no requirements on cell shapes. IsoAdvector is meant as 
a replacement for the interface compression with the MULES limiter implemented 
in the interFoam family of solvers.

The isoAdvector concept and code was developed at DHI and was funded by a Sapere
Aude postdoc grant to Johan Roenby from The Danish Council for Independent
Research | Technology and Production Sciences (Grant-ID: DFF - 1337-00118B - FTP).
Co-funding is also provided by the GTS grant to DHI from the Danish Agency for
Science, Technology and Innovation.

The ideas behind and performance of the isoAdvector scheme is documented in:

Roenby J, Bredmose H, Jasak H. 2016 A computational method for sharp interface 
advection. R. Soc. open sci. 3: 160405. 
[http://dx.doi.org/10.1098/rsos.160405](http://dx.doi.org/10.1098/rsos.160405)

Videos showing isoAdvector's performance with a number of standard test cases 
can be found in this youtube channel:

https://www.youtube.com/channel/UCt6Idpv4C8TTgz1iUX0prAA

## Requirements:

The isoAdvector code is developed and maintained in the newest OpenFOAM releases
but the script isoAdvector/bin/generateCodeForOldVersion copies the code in
isoAdvector/OpenFOAM to isoAdvector/OpenFOAM-[oldLoadedOFversion] and attempts 
to modify this code to become compatible with the older loaded OpenFOAM 
versions. 

A foam-extend version of the code is also available in isoAdvector/foam-extend. 
This version is not necessarily kept up to date.

## Installation:

0.  Source your OpenFOAM environment, e.g.: 

        source /home/$USER/OpenFOAM/OpenFOAM-4.1/etc/bashrc

1.  In a linux terminal download the package with git by typing:

        git clone https://github.com/isoAdvector/isoAdvector.git

2.  Execute the Allwmake script by typing (will finish in a ~1 min):

        cd isoAdvector
        ./Allwmake

    Applications will be compiled to your FOAM_USER_APPBIN and libraries will be
    compiled to your FOAM_USER_LIBBIN.
    
3.  Test installation with a simple test case by typing (finishes in secs):

        cp -r OpenFOAM-4.1/run/prescribedU/discInUniFlow/baseCase ~
        cd ~/baseCase
        ./Allrun
	
    Open Paraview and color the mesh by the alpha.water/alpha1 field.

4.  (Optional) If you want to test the code on unstructured meshes, a number of 
    such meshes can be donwloaded by using the downloadMeshes script in the bin 
    directory. These meshes will take up ~2GB and will be placed in a folder 
    called meshes. The meshes will be loaded by the scripts in discInUniFlows, 
    vortexSmearedDisc, sphereInUniFlow and smearedSphere cases in 
    run/isoAdvector.

    Alternatively, the meshes can be downloaded directly from here:

        http://dx.doi.org/10.5061/dryad.66840

    The downloaded meshes.tar.gz file should be extracted to the isoAdvector 
    root directory.
    
## Code structure:

`src/` 

* Contains the three classes implement the IsoAdvector advection scheme:
    - `isoCutFace`
    - `isoCutCell` 
    - `isoAdvection`
  These are compiled into a library named `libIsoAdvection`. 
* For comparison we also include the CICSAM, HRIC and mHRIC algebraic VOF 
  schemes in `finiteVolume` directory. These are compiled into a library called
  `libVOFInterpolationSchemes`.

`applications/` 

- `solvers/interFlow` 
    - A copy of the interFoam solver with the option to use isoAdvector instead
      of MULES in the interface advection step. To use isoAdvector add a 
      dictionary called isoAdvector to fvSolution with the contents:

      isoAdvector
      {
          //interfaceMethod can be set to "MULES" (default), "isoAdvector" or 
          //"fvSchemes". Use the latter option to use the HRIC, CICSAM or 
          //vofCompression schemes.

          interfaceMethod "isoAdvector";
          
          //isoFaceTol is the precision with wich the isosurface cutting a cell 
          //into subcells should reproduce the cell's volume fraction. Typically 
          //between 1e-6 and 1e-8 (default).

          isoFaceTol  1e-8;
          
          //surfCellTol defines which cells are treated as surface cells. If 
          //
          //  surfCellTol < alpha1 < 1 - surfCellTol
          //
          //a cell will be treated with isoAdvector. Typically between between 
          //1e-6 and 1e-8 (default).

          surfCellTol 1e-8;
          
          //nAlphaBounds is the number of times the volume preserving bounding 
          //procedure should be applied after the advection step to repair 
          //fluxes of unbounded cells. Default is 3.
          
          nAlphaBounds 3;

          //If snapAlphaTol > 0 then after advection and volume preserving 
          //bounding all remaining alpha's closer to 0 than snapAlphaTol will be 
          //set to 0 and all alpha's closer to 1 than snapAlphaTol will be set 
          //to 1. Default value is 1e-12.
          
          snapAlphaTol 1e-12; 

          //If clip is set to true/yes/1 then after advection and volume 
          //preserving bounding any alpha < 0 will be set to 0 and any alpha > 1 
          //will be set to 1.

          clip   true;

          //If prescribedU and PIMPLE.nCorrectors is set to -1, then the velocty
          //and pressure equations will not be solved. Useful for pure advection 
          //test cases.

          prescribedU true;

          //In cases with prescribed U there is an option to make the prescribed 
          //velocity field periodic by multiplying it by a factor 
          //cos(2*pi*runTime.time()/period) if period > 0:

          period      0;

          //In cases with prescribed U there is an option to reverse the 
          //velocity field when a when the time reverseTime is reached:

          reverseTime 0;
      }

      Please see cases in isoAdvector/run for examples of usage. Note that the 
      first versions of the isoAdvector code had the following other solvers:
      isoAdvector, mulesFoam, passiveAdvectionFoam. With the introduction of the 
      prescribedU switch described above, the functionallity of these solvers
      is now included in the interFlow solver and they have therefore been 
      removed.
- `utilities/preProcessing/isoSurf` 
    - Sets the initial volume fraction field for a sphere, a cylinder, a plane 
      or a sinus wave. See isoSurfDict for usage.
- `utilities/postProcessing/uniFlowErrors`
    - For cases with spheres and discs in steady uniform flow calculates errors 
      relative to exact VOF solution.
- `test/isoCutTester`
    - Application for testing isoCutFace and isoCutCell classes.

`run/`

- `prescribedU/` 
    - Contains pure advection test cases (prescribedU = true) from literature.
- `interFlow/` 
    - Contains test cases using interFlow coupling IsoAdvector with the PIMPLE 
      algorithm for the pressure-velocity coupling.
- `isoCutTester/`
    - Case for testing the isoCutFace and isoCutCell classes.
      
## Classes
	
###  IsoAdvection 

- Calculates the total volume of water, dVf, crossing each face in the mesh 
  during the time interval from time t to time t + dt.

### IsoCutCell

- Performs cutting of a cell given an isovalue and the alpha values interpolated 
  to the cell vertices.
- Calculates the submerged volume of a cell for a given isovalue.
- Calculates the isovalue given a specified target volume fraction of a cell.
  
### IsoCutFace

- Performs cutting of a face given an isovalue and alpha values interpolated to
  the face vertices.
- Calculates the submerged face area given an isovalue and alpha values 
  interpolated to the face vertices.

## Contributors:

* Johan Roenby <jro@dhigroup.com> (Inventor and main developer)
* Hrvoje Jasak (General coding guidance, consistent treatment of boundary faces 
  including processor boundaries, parallelisation, code clean up, provided 
  algebraic schemes, CICSAM, HRIC etc.)
* Henrik Bredmose (Participated in conceptual development)
* Vuko Vukcevic (Code review, profiling, porting to foam-extend, bug fixing, 
  testing)
* Tomislav Maric (Source file rearrangement)
* Andrew Heather (Code clean up, porting to OpenFOAM+)
