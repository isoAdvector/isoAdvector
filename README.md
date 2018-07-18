# Welcome to the IsoAdvector project

## What is IsoAdvector?

IsoAdvector is a geometric Volume-of-Fluid method for advection of a sharp 
interface between two incompressible fluids. It works on both structured and 
unstructured meshes with no requirements on cell shapes. IsoAdvector is meant as 
a replacement for the interface compression with the MULES limiter implemented 
in the interFoam family of OpenFOAM solvers.

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

## Compatibility:

The isoAdvector code is developed and maintained for the newest OpenFOAM 
releases but the script isoAdvector/bin/generateCodeForOldVersion copies the 
code in isoAdvector/OpenFOAM to isoAdvector/OpenFOAM-[oldLoadedOFversion] and 
attempts to modify this code to become compatible with an older sourced OpenFOAM 
version. The code should work with OpenFOAM-5.x, 4.1, 4.0, 3.0.1, 2.2.0,
v1606+ and v1612+. There are known to be issues with 2.2.2, 2.3.0, 2.3.1 and
2.4.0, so using these versions is not recommended.

OpenFOAM-v1706 and later vXXYY version contain an integration of the isoAdvector 
code with the interFlow solver named interIsoFoam. 
For more info on the v1706 integration see:

https://www.openfoam.com/releases/openfoam-v1706/numerics.php#numerics-isoadvector

Note that compiling the github.com/isoAdvector code with v1706 and later vXXYY 
versions may lead to problems because the OpenFOAM finiteVolume library already 
contains classes named isoAvection, isoCutCell and isoCutFace.

Compatibility with OpenFOAM-dev is regularly checked, but the developers tend to
frequently introduce minor API changes breaking the compatibility of isoAdvector
with the OpenFOAM-dev. If you experience compilation errors with dev, please 
notify me at the e-mail address below, and I'll fix it as soon as possible.

A foam-extend version of the code is also available in isoAdvector/foam-extend. 
This was developed for foam-extend-32 and will most likely need modifications to 
work with newer versions. It does not contain the latest code developments.

## Installation:

0.  Source your OpenFOAM environment, e.g.: 

        source /home/$USER/OpenFOAM/OpenFOAM-5.x/etc/bashrc

1.  In a linux terminal download the package with git by typing:

        git clone https://github.com/isoAdvector/isoAdvector.git

2.  Execute the Allwmake script by typing (will finish in ~1 minute):

        cd isoAdvector
        ./Allwmake

    Applications will be compiled to your FOAM_USER_APPBIN and libraries will be
    compiled to your FOAM_USER_LIBBIN (You can check where e.g. the environmental 
    variable FOAM_USER_APPBIN point to by typing 

        echo $FOAM_USER_APPBIN

    in the terminal).
    
3.  Test installation with a simple test case by typing (finishes in secs):

        cp -r OpenFOAM-5.x/run/prescribedU/discInConstantFlow/baseCase ~
        cd ~/baseCase
        ./Allrun
	
    Open Paraview and color the mesh by the alpha.water/alpha1 field.

4.  (Optional) If you want to test the code on unstructured meshes, a number of 
    such meshes can be donwloaded by using the downloadMeshes script in the bin 
    directory. These meshes will take up ~2GB and will be placed in a folder 
    called meshes. The meshes will be loaded by the scripts in 
    discInConstantFlows, discInReversedVortexFlow, sphereInConstantFlow and 
    sphereInReversedVortexFlow cases in run/prescribedU.

    Alternatively, the meshes can be downloaded directly from here:

        http://dx.doi.org/10.5061/dryad.66840

    The downloaded meshes.tar.gz file should be extracted to the isoAdvector 
    root directory.

    Unstrictured meshes can also be generated with the scripts generateTetUnitCube
    and generateTriUnitSquare in the bin folder. This requires gmsh (see headers).
    
## Code structure:

`src/` 

* Contains the three classes implement the IsoAdvector advection scheme:
    - `isoCutFace`
    - `isoCutCell` 
    - `isoAdvection`
  These are compiled into a library named `libIsoAdvection`. 
* For comparison we also include the CICSAM, HRIC and mHRIC algebraic VOF 
  schemes in `finiteVolume` directory. These were previously compiled into a 
  library called `libVOFInterpolationSchemes` but are currently not compiled
  and are only included as legacy code.

`applications/` 

- `solvers/interFlow` 
    - A copy of the interFoam solver with the option to use isoAdvector instead
      of MULES in the interface advection step. To use add set application in 
      controlDict to interFlow and add the following to the dicitonary
      fvSolution.solvers."alpha.water.*":

      ```
      "alpha.water.*"
      {
          //interfaceMethod can be set to "MULES" (default), or "isoAdvector".

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

          //If snapTol > 0 then after advection and volume preserving 
          //bounding all remaining alpha's closer to 0 than snapTol will be 
          //set to 0 and all alpha's closer to 1 than snapTol will be set 
          //to 1. Default value is 1e-12.
          
          snapTol 1e-12; //Previously called clipAlphaTol (and snapAlphaTol)

          //If clip is set to true/yes/1 then after advection and volume 
          //preserving bounding any alpha < 0 will be set to 0 and any alpha > 1 
          //will be set to 1.

          clip   true; //Previously called clipAlpha

          //To write out case/isoFaces/isoFaces#tIndex.obj change this to true:

          writeIsoFaces false;

          //For tri and tet meshes the standard isoAdvector method may result in 
          //large variations in the interface normal orientation in neighbouring
          //cells. A much smoother interface normal orientation is obtained by 
          //enforcing use of a smoothed gradient for the isoface orientations.
          //This method is activated by changing this to true:

          gradAlphaNormal false;
      }
      ```

      To use the functionality of a prescribed flow you can set 
      momentumCorrector to no and nCorrectors to -1 in the fvSolution.PIMPLE 
      directory and add the following to the dicitonary fvSolution.solvers.U:

      ```
      U
      {
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
      ```

          
      Please see cases in OpenFOAM/run for examples of usage. Note that the 
      first versions of the isoAdvector code had the following other solvers:
      isoAdvector, mulesFoam, passiveAdvectionFoam. With the introduction of the 
      prescribedU switch described above, the functionallity of these solvers
      is now included in the interFlow solver and they have therefore been 
      removed.
- `utilities/preProcessing/setAlphaField` 
    - Sets the initial volume fraction field for a sphere, a cylinder, a plane 
      or a sinus wave. See setAlphaFieldDict for usage. Previsouly called 
      isoSurf.
- `utilities/postProcessing/calcAdvectErrors`
    - For cases with spheres and discs in steady uniform flow calculates errors 
      relative to exact VOF solution. Previously called uniFlowErrors.
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

* Johan Roenby, STROMNING, <johan@stromning.com> (Inventor and main developer)
* Hrvoje Jasak, University of Zagreb (General coding guidance, consistent treatment 
  of boundary faces including processor boundaries, parallelisation, code clean up, 
  provided algebraic schemes, CICSAM, HRIC etc.)
* Henrik Bredmose, DTU Wind Energy (Participated in conceptual development)
* Vuko Vukcevic, University of Zagreb (Code review, profiling, porting to 
  foam-extend, bug fixing, testing)
* Andrew Heather, OpenCFD (Code clean up, porting to OpenFOAM+)
* Henning Scheufler, DLR (Extensive validation, gmsh based mesh generation)
