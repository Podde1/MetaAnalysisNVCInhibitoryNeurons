function [lb ,ub] = setLBandUB(thetaOptParameterNames)
%% Generate parameter bounds
lb = -5*ones(size(thetaOptParameterNames));
ub = 5*ones(size(thetaOptParameterNames));

% thetaOptParameterNames = thetaOptParameterNames(OptParamBool); % Names of parameters that are included in parameter estiamtion.

%% set specific parameter estimation bounds for specific paramters (where applicable)
ub(ismember(thetaOptParameterNames, {'k_u1' 'k_u2' 'k_u3' 'ku4'})) = 5; 
lb(ismember(thetaOptParameterNames, {'k_u1' 'k_u2' 'k_u3' 'ku4'})) = -5; 
% ub(ismember(thetaOptParameterNames, {'N_NObasal' 'N_NPYbasal' 'N_Pyrbasal'})) = 3; 
% lb(ismember(thetaOptParameterNames, {'N_NObasal' 'N_NPYbasal' 'N_Pyrbasal'})) = -3; 
ub(ismember(thetaOptParameterNames, {'N_NOMax' 'N_NPYMax' 'N_PyrMax' 'N_SOMMax'})) = 5; 
lb(ismember(thetaOptParameterNames, {'N_NOMax' 'N_NPYMax' 'N_PyrMax' 'N_SOMMax'})) = -3; 
% lb(ismember(thetaOptParameterNames, {'N_SOMMax'})) = 0; 

% lb(ismember(thetaOptParameterNames, {'N_NPYMax' 'N_PyrMax'})) = 0; 

% ub(ismember(thetaOptParameterNames, {'kPF1' 'kPF2' 'kIN' 'kIN2' 'kINF' 'kINF2'})) = 5;
% lb(ismember(thetaOptParameterNames, {'kPF1' 'kPF2' 'kIN' 'kIN2' 'kINF' 'kINF2'})) = -5;
% 
% ub(ismember(thetaOptParameterNames, {'kPF1' 'kPF2' 'kPF3' 'kPF4' 'kPF5' 'kIN' 'kIN2' 'kIN3' 'kIN4' 'kIN5'})) = 5;
% lb(ismember(thetaOptParameterNames, {'kPF1' 'kPF2' 'kPF3' 'kPF4' 'kPF5' 'kIN' 'kIN2' 'kIN3' 'kIN4' 'kIN5'})) = -5;

lb(ismember(thetaOptParameterNames, {'Km' 'Km2'})) = -12;  % Km saturation parameters, allowed to be small
lb(ismember(thetaOptParameterNames, {'sinkN_NPY' 'sinkN_Pyr' 'sinkN_Astro' 'sinkN_SOM' 'sinkCa_NO' 'sinkCa_NPY' 'sinkCa_Pyr' 'sinkCa_Astro' 'sinkCa_SOM'})) = log10(1/0.78);    %ksink for N and Ca2+	  
lb(ismember(thetaOptParameterNames, {'sinkN_Astro'})) = -5; % sinkN_Astro
lb(ismember(thetaOptParameterNames, {'sinkN_NO'})) = 0;     % sinkNO 
lb(ismember(thetaOptParameterNames, {'sinkN_SOM'})) = -5;   % sinkN_SOM

% lb(ismember(thetaOptParameterNames, {'nNO'})) = log10(0.1);
% ub(ismember(thetaOptParameterNames, {'nNO'})) = log10(4);

lb(ismember(thetaOptParameterNames, {'kCMRO2_Pyr' 'kCMRO2_SOM'})) = -2;
lb(ismember(thetaOptParameterNames, {'kCMRO2_NO' 'kCMRO2_NPY','kCMRO2_Astro'})) = -5;

% Circuit parameters
lb(ismember(thetaOptParameterNames, {'K1'})) = 0.001;     
lb(ismember(thetaOptParameterNames, {'K2'})) = log(1.1);          
lb(ismember(thetaOptParameterNames, {'K3'})) = 0;       % K3 continuously ends up at lower bound = 2 (default). Hence lb set lower         
ub(ismember(thetaOptParameterNames, {'K1'})) = log(2);     
ub(ismember(thetaOptParameterNames, {'K2'})) = log(2);          
ub(ismember(thetaOptParameterNames, {'K3'})) = 8;

lb(ismember(thetaOptParameterNames, {'vis1'})) = 0;     
lb(ismember(thetaOptParameterNames, {'vis2'})) = 0; % vis2 continuously ends up at lower bound = 1 (default). Hence lb set lower         
lb(ismember(thetaOptParameterNames, {'vis3'})) = 1;
ub(ismember(thetaOptParameterNames, {'vis1'})) = 2;     
ub(ismember(thetaOptParameterNames, {'vis2'})) = 3;          
ub(ismember(thetaOptParameterNames, {'vis3'})) = 5; % vis3 continuously ends up at upper bound = 3 (default). Hence ub set higher

lb(ismember(thetaOptParameterNames, {'k_capperm'})) = log10(0.1);
ub(ismember(thetaOptParameterNames, {'k_capperm'})) = log10(10); 

lb(ismember(thetaOptParameterNames, {'k_leak'})) = -3;
ub(ismember(thetaOptParameterNames, {'k_leak'})) =  2; 

lb(ismember(thetaOptParameterNames, {'k_yBOLD', 'k_yBOLD2'})) = -3;

end