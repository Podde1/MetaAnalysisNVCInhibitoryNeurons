function [] = GenerateModels()
%% check checkCompatibility 
if ~contains(version,'R2017b')
    warning('The following script require MATLAB release R2017b! Current version is: %s', version)
end 
if ~contains(path,'AMICI')
    warning('"AMICI" is not found in matlab''s search path! The AMICI-toolbox might not be installed.')
end 
   

%% compile models
amiwrap('Sten2023_v21','Sten2023_v21_syms','.\Models\mex_and_simFiles\')
amiwrap('Sten2023_v21_steadyState','Sten2023_v21_steadyState_syms','.\Models\mex_and_simFiles\')

amiwrap('Sten2023_v31','Sten2023_v31_syms','.\Models\mex_and_simFiles\')
amiwrap('Sten2023_v31_steadyState','Sten2023_v31_steadyState_syms','.\Models\mex_and_simFiles\')

amiwrap('Sten2023_v21_BOLDcontributions','Sten2023_v21_BOLDcontributions_syms','.\Models\mex_and_simFiles\')
amiwrap('Sten2023_v21_BOLDcontributions_steadyState','Sten2023_v21_BOLDcontributions_steadyState_syms','.\Models\mex_and_simFiles\')

disp('---------- ALL MODELS COMPLIED SUCCESSFULLY ----------')
end

