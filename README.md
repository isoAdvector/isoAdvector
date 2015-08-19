# IsoAdvection Project  

This is the repository of the IsoAdvection project. 

Developed by Johan Rønby <jro@dhigroup.com>  

Contributors:  

Tomislav Maric <tomislav@sourceflux.de> 

# Possible future issues 

Add feature and improvement ideas/requests here. Some of them will be taken over by the development team and converted to issues. 

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
