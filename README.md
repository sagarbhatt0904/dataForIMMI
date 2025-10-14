# Data for the paper 
## *"Modeling the Impact of Grain Morphology and Porosity on Creep Anisotropy in Laser Powder Bed Fusion Materials"*

This repository contains the data, input files, and results for the simulations presented in the paper. These simulations were conducted using the MOOSE based CPFEM framework DEER (https://github.com/Argonne-National-Laboratory/deer), which utilizes the constitutive material models from NEML (https://github.com/Argonne-National-Laboratory/neml)

- `common/` contains 
    - the  `<load-dir>.i` MOOSE input files for DEER 
    - the `.xml` model file for NEML
    - an example script to run the simulations
- `microstructureData/` contains
    - the simulated microstructure data in a `.txt` file
    - a `.csv` file where each line contains the Bunge euler angles for the grain ID matching the line number
    - a `.e` exodus mesh file of generated from the microstructure
    - a `.tex` orientation file containing orientations for each grain in the mesh
    - a python script to generate the `.tex` file for the mesh from the microstructure data
- `simulationResults/` contains the creep results. 
    - the top level folder names indicate the loading directions
    - under each top level folder, the subfolder names indicate the initial cavity size $a_0$
    - in each subfolder, there is a `.csv` file that contains the results from the simulations and a `gbProp.i` file that contains the parameters for the grain boundary cavitation model.
- `tensileTests` contains two sets of simulations in their respective folders -- `withCavitation` (with grain boundary cavitation model turned on) and `withoutCavitation` (no grain boundary cavitation). Each folder contains:
    - `common` - directory with `<load_dir>.i` MOOSE input files for DEER and `.xml` model file for NEML
    - `load_dir/` - directory with a `.csv` file containing results from tensile simulations.
    - under `withCavitation`, for the loading direction `z`, the input files and result directory names are appended with the cavity size used to study sensitivity of cavity sizes on yield strength.

