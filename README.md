# Energy-verification-for-OpenHFDIB-DEM

This repository contains a series of test cases with a verification tool developed for evaluating energy conservation in CFD-DEM simulations using open-source OpenHFDIB-DEM (https://github.com/techMathGroup/openHFDIB-DEM) project.

The project is based on OpenFOAMv8 (https://openfoam.org), combining Computational Fluid Dynamics (CFD) and the Discrete Element Method (DEM) with a Hybrid Fictitious Domain-Immersed Boundary (HFDIB) approach.

In pure DEM, results from openHFDIB-DEM can be compared with simulations performed using the open-source DEM solver LIGGGHTS (https://github.com/CFDEMproject/LIGGGHTS-PUBLIC).

Code capabilities
-----------------
* Verifies energy conservation in pure DEM and CFD-DEM simulations.
* Tracks various energy components, such as kinetic, potential, rotational, dissipated, contact, drag and fluid energy. 
* Supports both spherical and arbitrary shaped particle cases.

Usage
-----------------
For simulation setup and execution, please refer to the official OpenHFDIB-DEM repository.

Once a simulation has completed, run verification using:

* python3 evaluateSimulationCompleteData.py

This script will process the simulation log file, compute all energy components, and generate postProc-results/ folder containing files with energy data. 