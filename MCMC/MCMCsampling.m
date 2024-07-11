function [] = MCMCsampling(theta,modelName,ExperimentIdx,variables2Fit)
% theta = load(Parameters/...)
% modelName = 'Sten2023_v21'; 
% ExperimentIdx = [1 4 5 8:15]; % specify the experiment index for which experiment to run parameter estimation for 
% variables2Fit = {'CBF','BOLD','HbT','HbO','HbR'}; % specify the name as a string of which varaibles i the experiment to fit o during parameter estimation

addpath(genpath('./lib'))
addpath(genpath('./Data'))
addpath(genpath('./Models'))

% run model set up function 
[circuitConstants, DataStructure, model, variableNames, saveString, modelObservablesIndex, parameterMappingTable, options] = modelSetup(modelName,ExperimentIdx,variables2Fit);

%% Create file where all chi-2 acceptable parameters are saved
FileName=sprintf('MCMC_%s_%s.dat',modelName,datestr(now,'yymmdd_HHMMSS')); %use for real run.
% FileName = 'TempTESTFile.dat'; % use for testing and debugging such that a new file is not created for each run.
FID = fopen(FileName,'wt');


%% Set parameter bounds 
thetaOptParameterNames = model.parameters;

%% Ensure that the theta is as long as the model parameter vector. pad with zeros if that is not the case
if numel(model.parameters)>numel(theta) % if the model parameter vector is longer than the loaded theta 
    nParameterDiff = numel(model.parameters)-numel(theta); % calulate the difference in parameter vector length 
    theta(end+1:end+nParameterDiff) = zeros(1,nParameterDiff); % add zeros to the unkown parameter values 
end

%% configure the total parameter estimation vector to contian different parameter sets for different experiments 
for i = 1:size(parameterMappingTable,1) %loop throught the parameter Mapping table.
    if any(parameterMappingTable.TotalParameterIndex{i}>numel(theta)) % if any of the TotalParameterIndex on th i:th row exceed the number of parameters in the loaded theta-vector. 
        theta(parameterMappingTable.TotalParameterIndex{i}) = theta(parameterMappingTable.ModelParameterIndex{i}); %set the inital value of any not previously existing parameter values to that of the corresponding parameter index .
        thetaOptParameterNames(parameterMappingTable.TotalParameterIndex{i}) = thetaOptParameterNames(parameterMappingTable.ModelParameterIndex{i}); % set the corresponding parameters names for the added parameters.
    end

    if any(parameterMappingTable.TotalParameterIndex{i}>numel(thetaOptParameterNames)) % if any of the TotalParameterIndex on th i:th row exceed the number of parameters to optimize 
        thetaOptParameterNames(parameterMappingTable.TotalParameterIndex{i}) = thetaOptParameterNames(parameterMappingTable.ModelParameterIndex{i}); % set the corresponding parameters names for the added parameters.
    end
end

%% Allow the user to specify the index of specific parameters to estimate. All other parameters will be fixated at the inital value.  
OptParamBool = false(1, numel(thetaOptParameterNames)); %preallocate a boolean value for every parameter. The boolean value is used to determine if a parameter is included in the parameter estimation or if it is keep constant. 1 = included, 0 = keept constant. 

maxNParametersInSetup = max(max(cell2mat(parameterMappingTable.TotalParameterIndex.')),numel(model.parameters)); %compare and get the larger numebr of the largest index in the totalParameter Idx or the number of parameters in the model.
OptParamIndexes = 1:maxNParametersInSetup; %set index of the parameters to be included in parameter estimation. Default is the full maximal number of parameters in the current setup , depending on which experiments are choosen. 

if strcmp("Sten2023_v31", modelName)
    OptParamIndexes = [1:80 83 86 87 90];
else
    OptParamIndexes = [];
end

OptParamBool(OptParamIndexes) = true; %set parameters to optimise to have a true flag

thetaOpt = theta(OptParamBool).'; % build the parameter vector that is passed to parameter estimation algorithm.
thetaConstants = theta(~OptParamBool).';% build the parameter vector that is keept constant during parmater estimation.

%% Generate parameter bounds
thetaOptParameterNames = thetaOptParameterNames(OptParamBool); % Names of parameters that are included in parameter estiamtion.
[lb,ub] = setLBandUB(thetaOptParameterNames);

parameters.min = lb.';
parameters.max = ub.';

parameters.name = thetaOptParameterNames;
parameters.number = length(thetaOpt);           % number of parameters

%% Set the objective function
proxyObjectiveFunction = @(X) objectiveFunction(X, modelName, thetaConstants, OptParamBool, circuitConstants, DataStructure,ExperimentIdx,variableNames,modelObservablesIndex, parameterMappingTable, model, options,0,FID);

%% Options
optionsPesto = PestoOptions(); % loads optimization options and options for everything basically

optionsPesto.obj_type = 'negative log-posterior';  % Works with least squares cost which we typically been using, i.e the goal is the minimize. so keep it as is
optionsPesto.n_starts = 1; % number of multi-starts optimization performed <-- ignore this, better to start sampling in a previously estimated optima
optionsPesto.mode = 'visual'; %'silent'' it will plot the results itteratively. 

optionsPesto.comp_type = 'sequential'; % no parallel computing 
%% If using SBToolbox, it won't return everything as needed (i.e hessian and derivatives, have this line active)
% The algorithm does not need any sensitivities to work, therefore only one output
optionsPesto.objOutNumber=1;

%% Markov Chain Monte Carlo sampling -- Parameters
% Values for the parameters are sampled by using an advanced Parallel Tempering (PT)
% algorithm. This way, the underlying probability density of the parameter 
% distribution can be captured. if only ONE temperature is used, this is
% effectively an adapted Metropolis algorithm single-chain algorithm.

% Building a struct covering all sampling options:
optionsPesto.MCMC = PestoSamplingOptions();
optionsPesto.MCMC.nIterations = 1e5; 
optionsPesto.MCMC.mode = optionsPesto.mode;
%% Using RAMPART
optionsPesto.MCMC.samplingAlgorithm     = 'RAMPART';
optionsPesto.MCMC.RAMPART.nTemps           = size(thetaOpt,1);
optionsPesto.MCMC.RAMPART.exponentT        = 1000;
optionsPesto.MCMC.RAMPART.maxT             = 2000;
optionsPesto.MCMC.RAMPART.alpha            = 0.51;
optionsPesto.MCMC.RAMPART.temperatureNu    = 1e3;
optionsPesto.MCMC.RAMPART.memoryLength     = 1;
optionsPesto.MCMC.RAMPART.regFactor        = 1e-8;
optionsPesto.MCMC.RAMPART.temperatureEta   = 10;

optionsPesto.MCMC.RAMPART.trainPhaseFrac   = 0.2;
optionsPesto.MCMC.RAMPART.nTrainReplicates = 5;

optionsPesto.MCMC.RAMPART.RPOpt.rng                  = 1;
optionsPesto.MCMC.RAMPART.RPOpt.nSample              = floor(optionsPesto.MCMC.nIterations*optionsPesto.MCMC.RAMPART.trainPhaseFrac)-1;
optionsPesto.MCMC.RAMPART.RPOpt.crossValFraction     = 0.2;
optionsPesto.MCMC.RAMPART.RPOpt.modeNumberCandidates = 1:20;
optionsPesto.MCMC.RAMPART.RPOpt.displayMode          = 'text';
optionsPesto.MCMC.RAMPART.RPOpt.maxEMiterations      = 100;
optionsPesto.MCMC.RAMPART.RPOpt.nDim                 = parameters.number;
optionsPesto.MCMC.RAMPART.RPOpt.nSubsetSize          = 1000;
optionsPesto.MCMC.RAMPART.RPOpt.lowerBound           = parameters.min;
optionsPesto.MCMC.RAMPART.RPOpt.upperBound           = parameters.max;
optionsPesto.MCMC.RAMPART.RPOpt.tolMu                = 1e-4 * (parameters.max(1)-parameters.min(1));
optionsPesto.MCMC.RAMPART.RPOpt.tolSigma             = 1e-2 * (parameters.max(1)-parameters.min(1));
optionsPesto.MCMC.RAMPART.RPOpt.dimensionsToPlot     = [1,2];
optionsPesto.MCMC.RAMPART.RPOpt.isInformative        = [1,1,ones(1,optionsPesto.MCMC.RAMPART.RPOpt.nDim-2)];
      
optionsPesto.MCMC.theta0 = thetaOpt;
optionsPesto.MCMC.sigma0 = optionsPesto.MCMC.nIterations * eye(numel(thetaOpt));

%% Run the sampling
warning('off','all')
parameters = getParameterSamples(parameters, proxyObjectiveFunction, optionsPesto);
warning('on','all')

cd('../../')
%% Create folder with name MCMC_DAY-MONTH-YEAR hour-min-sec (Very advanced matlab)

folderStr = [datestr(now,'yymmdd_HHMMSS'),'_',modelName];
folderStr = strrep(folderStr,{':',' '},'_');
folderStr = folderStr{2};
mkdir(['./Results/MCMCResult_' folderStr])
FolderName=fullfile(['./Results/MCMCResult_' folderStr]);
%% Save results to folder
save(fullfile(FolderName,'parameters.mat'),'parameters','theta','parameterMappingTable','saveString')
fclose(FID);
movefile(sprintf('./Models/%s/%s',modelName,FileName),FolderName);
%% Takes all open figures and saves them
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    FigName = get(FigHandle, 'Name');
    try
        savefig(FigHandle, fullfile(FolderName,['fig', num2str(iFig), '.fig']));
    catch
    end
end

% extract max/min results from MCMC
ExtractMCMC(fullfile(FolderName,'parameters.mat'))
