# Wave Scattering & Geometrical Optics Simulation
This project analyzes electromagnetic wave behavior using both numerical and analytical techniques.  
It includes two main components:

1. **Electromagnetic scattering from two dielectric cylinders** using the Volume Method of Moments (VMoM), with comparison to Born and Rytov approximations.

2. **Ray propagation in an inhomogeneous medium** using Geometrical Optics (GO), including ray trajectories, phase accumulation, turning points, caustics, and 2D field mapping.

All simulations are implemented in MATLAB and demonstrate the combination of EM theory, numerical methods, and visualization tools.

##  How to Run the Code

1. Open **MATLAB**.
2. Add the project folder to your MATLAB path:
   ```matlab
   addpath(genpath('path/to/project'))
3.Run the main script for the desired part:
##  Project Structure-
Par 1 :scattering simulation

1.**section 1:
1.1. **Volume Method of Moments implementation
sec1_1.m
1.2. **Born and Rytov Approximations 
 sec1_2.m
2. **section 2: Inverse Scattering Using Born Approximation
 ****sec_2.m
Part 2: Geometrical Optics implementation
│ ├── part2_RayPaths.m
│ ├── part2_phase.m
│ ├── part2_caustics.m
│ └──part2_Electric_Field.m






