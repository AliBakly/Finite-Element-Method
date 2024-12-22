# Thermo-Mechanical Gripper FEM Project

This repository contains a Finite Element Method (FEM) project for simulating **thermal and mechanical** behavior of a gripper device. The project is associated with the course **FHLF25, Finite Element Method and Introduction to Strength of Materials** at Lund University (LTH). 

![Thermo-Elastic Deformation](ThermoDeformation.gif)

## Contents

- **`FEMProject.m`**  
  Main MATLAB script implementing the following tasks:
  1. **Stationary Heat Convection**: Solves for the steady-state temperature distribution with convection and heat flux boundary conditions.
  2. **Transient Heat Convection**: Time-dependent thermal analysis (both snapshots and an animation).
  3. **Thermoelasticity (Transient)**: Uses the time-evolving temperature fields to solve for mechanical displacements and stresses (von Mises).
  4. **Visualization**: Generates plots of temperature distributions, von Mises stress fields, and a short “movie” of the deformation over time.
  
- **`FEM.pdf`**  
  A more detailed **report** describing the theoretical background, problem setup, derivations, and results.

- **`4nmesh.mat`** / **`highermesh.mat`**  
  Example meshes exported from MATLAB’s PDETool.  
  - `4nmesh.mat` is a coarser mesh of the geometry.  
  - `highermesh.mat` is a finer mesh used to compare results.

- **`calfem/`**  
  Directory containing or referencing the [CALFEM](https://github.com/CALFEM/calfem-matlab) library code. CALFEM is used for element routines such as `flw2te`, `plante`, `plantml`, etc.

- **`plantml.m`**  
  A local function (if not inside `calfem/`) used to compute element capacity matrices for transient heat conduction.  

- **`README.md`**  
  This file. Provides instructions, references, and context for the project.

## Requirements

1. **MATLAB** (or compatible environment such as Octave, though not tested).  
2. **CALFEM for MATLAB**  
   - This project relies on CALFEM routines (`flw2te`, `plante`, `plantml`, etc.).  
   - You can download CALFEM from [CALFEM/matlab on GitHub](https://github.com/CALFEM/calfem-matlab).  
   - Ensure that the path to CALFEM is added via `addpath(genpath("calfem/fem/"))`, or adjust accordingly in the code.

3. **Mesh Files**  
   - `4nmesh.mat` or `highermesh.mat` which contain geometry and boundary labeling from PDETool.

## How to Run

1. **Clone or download** this repository.
2. **Install CALFEM** (if you have not done so already), and confirm you can import it in MATLAB.
3. **Open `FEMProject.m`** in MATLAB.
4. Press **Run** (or call `FEMProject` from the command line). 

   - The script will load `4nmesh.mat` by default.  
   - It will build global stiffness and capacity matrices, solve the conduction problem, and then proceed with the thermoelastic calculations.
   - Finally, it will produce multiple **figures** for temperature fields, von Mises stress, and a short “movie” of the deforming gripper under thermal expansion.

5. **Check the Output**  
   - You should see figures with temperature distributions (stationary and transient), snapshots in time, an animated conduction “video,” and a final displacement + von Mises plot.

## Project Outline

1. **Stationary Heat Convection (Task a)**
   - Builds a global conduction matrix **K** and force vector **f** from conduction, convection, and flux boundary conditions.  
   - Solves the linear system **K * T = f** to obtain nodal temperatures.

2. **Transient Heat Convection (Task b)**
   - Constructs a global capacity matrix **C** for the time-stepping scheme.
   - Advances the temperature field in time by solving:
     **(C + Δt*K) * T(n+1) = Δt*f + C*T(n)**  
   - Produces snapshots at select time intervals, plus an animation of temperature evolution.

3. **Thermoelasticity (Task c)**
   - Uses the temperature field from Task b to compute thermal strains.
   - Assembles plane-strain elasticity equations:
     **K_el * u = f_thermal**  
   - Solves for displacements u at each time step, then calculates the final von Mises stress distribution.

4. **Visualization**
   - Plots the final stress distribution (highlighting the node with maximum stress).
   - Animates the combined temperature and displacement fields over time.
  

## Course and Credits

- This project is part of **FHLF25 – Finite Element Method and Introduction to Strength of Materials** at **Lund University (LTH)**.
- In collaboration with **Davy Than**.
- The PDETool geometry/mesh and boundary definitions are from the gripper geometry pictured in the report.

## License / Citation

- If you use or extend this code, please credit the authors and the [CALFEM library](https://github.com/CALFEM/calfem-matlab).
- The project is for educational purposes within **FHLF25**.
