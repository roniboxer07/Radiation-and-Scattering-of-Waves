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
3.Run the main script for the desired part

## Project Structure

### Part 1: Electromagnetic Scattering Simulation

**1. Section 1**  
&nbsp;&nbsp;&nbsp;**1.1. Volume Method of Moments implementation**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;• `sec1_1.m`  
&nbsp;&nbsp;&nbsp;**1.2. Born and Rytov Approximations**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;• `sec1_2.m`

**2. Section 2: Inverse Scattering Using Born Approximation**  
&nbsp;&nbsp;&nbsp;• `sec_2.m`

---

### Part 2: Geometrical Optics Implementation

&nbsp;&nbsp;&nbsp;**1.1.** `part2_RayPaths.m`  
&nbsp;&nbsp;&nbsp;**1.2.** `part2_phase.m`  
&nbsp;&nbsp;&nbsp;**1.3.** `part2_caustics.m`  
&nbsp;&nbsp;&nbsp;**1.4.** `part2_Electric_Field.m`

## Requirements
MATLAB R2021a or newer, Basic MATLAB linear algebra support, MATLAB Graphics toolbox.





