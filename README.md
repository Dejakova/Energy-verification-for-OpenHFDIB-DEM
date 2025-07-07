# Energy-verification-for-OpenHFDIB-DEM

This repository contains a series of test cases with a verification tool developed for evaluating energy conservation in CFD-DEM simulations using open-source OpenHFDIB-DEM (https://github.com/techMathGroup/openHFDIB-DEM) solver.

The solver is based on OpenFOAM (https://openfoam.org), combining Computational Fluid Dynamics (CFD) and the Discrete Element Method (DEM) with a Hybrid Fictitious Domain-Immersed Boundary (HFDIB) approach.

Code capabilities
-----------------
* Verifies energy conservation in pure DEM and CFD-DEM simulations.
* Tracks various energy components, such as kinetic, potential, rotational, dissipated, contact, drag and fluid energy. 
* Supports both spherical and arbitrary shaped particle cases.

Usage
-----------------
For simulation setup and execution, please refer to the official OpenHFDIB-DEM repository.

Once a simulation has completed, you run verification using:

* python3 evaluateSimulationCompleteData.py

This script will process the simulation log file, compute all energy components, and generate postProc-results/ folder containing files with energy data. 