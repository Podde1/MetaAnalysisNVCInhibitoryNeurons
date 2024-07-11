## Prerequisites
To run this repo, the following is required:

Matlab R2017b(used to compile the model files), with add ons:
    - mex compiler (preferably MinGw-64): https://www.mathworks.com/matlabcentral/fileexchange/52848-matlab-support-for-mingw-w64-c-c-compiler
    - Symbolic Math Toolbox: https://www.mathworks.com/products/symbolic.html

For full code comparability Matlab >=2023a (to run scripts)

Advanced Multilanguage Interface to CVODES and IDAS https://github.com/AMICI-dev/AMICI. Under the "lib" folder a slightly altered version is available, where issues related to the compilation of events are fixed. The mex compilation of models can cause errors if the path is to long, if that is the case, try moving the AMICI folder further up in the directory tree.

PESTO: Parameter EStimation TOolbox https://github.com/ICB-DCM/PESTO

Furthermore, to run the optimization the MEIGO toolbox in needed: https://bitbucket.org/jrbanga_/meigo64/src/master/

## Get the results from the paper
Simply run the file `plotAllResults.m` to get all figures from the paper. 

## Parameter estimation
To perform parameter estimations, the "parameterEstimation/" folder has prepared scripts. 
`runParallelParameterEstimation.m` requires the `Parallel Computing Toolbox` (matlab add on)
`parameterEstimation` run the parameter estimation in sequence instead of parallel and do not require the add on

## Uncertainty collection 
Prediction Profile Likelihood (PPL) and reverse PL (revPL) can be run from the "ProfileLikelihoodAnalysis/" Folder 

## Structure of the repository
Data/ - the experimental data is stored here

Figures/ - preallocated storage for figures to be saved here

lib/ - contains various functions and tools for the project

lib/ObjectiveFunction - the objective functions used for the various parameter estimations

lib/plotFunctions - all plot scripts for the project

lib/simulationFunctions - the simulation function for both models

Models/symsFiles - the model files (_syms) is stored here

Models/mex_and_simFiles - after the model is compiled, the .mex file and Simulate.m file is stored here and is called by the simulations function

Parameters/ - All prepared results to plot the published results (figures) are stored here

Results/ - preallocated storage of the various results if parameter estimation/uncertainty is performed