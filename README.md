# Overview
This repository contains the source code, data files and parameters for the paper ``A Novel Decomposition-Based Exact Solution Approach for Lot-sizing and Scheduling Problem" by Eyüp Ensar Işık, Z. Caner Taşkın and Semra Ağralı.

# Dependencies
* Microsoft Visual Studio 2022 with C++ language support
* IBM ILOG CPLEX 22.1 Windows x64 version
# File structure
* Data/: Data files. Each problem instance has a corresponding .dat file. Benchmark instances taken from James and Almada-Lobo (2011). 
* Code: Source code 
* Code\VS2022: Visual Studio project and solution files
* Code\Run: Batch files to run experiments for the Decomposition Algorithm and MIP models 
* Code\Parameters: Parameter files for computational experiments

# Build and run instructions
* Clone repository to a local folder. Ensure that dependencies are installed in their default locations (if not, adjust Visual Studio project settings to match the installed location and CPLEX version).
* Open Code\VS2022\GLSP.sln file in Visual Studio.
* Choose ``Release" build configuration and build the solution.
* Navigate to Code\Run on command line
* Run singlemachineStd.bat, singlemachineNF.bat, singlemachineCC.bat, singlemachineDA.bat to run computational experiments reported in the paper. 
* Outputs will be generated in Code\VS2022\Results folder. 

# References
* James RJ and Almada-Lobo B. Single and parallel machine capacitated lotsizing and scheduling: New iterative MIP-based neighborhood search heuristics. Computers & Operations Research 2011;38:1816–25.
