function [] = PPL(theta,resultFolderName, modelName,ExperimentIdx,variables2Fit,Index,referenceValue,step,polarity)
% modelName - string with the specific model name "_syms", "mex*" appendix not included.
% ExperimentIdx - specify the experiment index for whihc experiment to run parameter estimation for 
% variables2Fit - specify the name as a string of which varaibles i the experiment to fit o during parameter estimation
% Index - parameterindexes to evaluate 
% referenceValue - value to set and have as reference
% step - stepsize
% Polarity - direction (1, -1) to evaluate towards

addpath(genpath('./lib'))
addpath(genpath('./Data'))
addpath(genpath('./Models'))

%% Set up a result folder 
if isempty(resultFolderName) % if resultFolderName is empty create a template folder
    resultFolderName = char(strcat(modelName, '_Experiment_',strjoin(string(ExperimentIdx),'_'))); % create a result folder name based on hte experimetIndex and the model name. 
end

ResultLocation=['./Results/' resultFolderName,'/'];

if ~exist(ResultLocation,'dir') % if the location does not exist create it 
    mkdir(ResultLocation)
end
%% Run inital model set up
[circuitConstants, DataStructure, model, variableNames, saveString, modelObservablesIndex, parameterMappingTable, options] = modelSetup(modelName,ExperimentIdx,variables2Fit);

thetaOptParameterNames = model.parameters;
%% Ensure that the theta is as long as the model parameter vector. pad with zeros if that is not the case
if numel(model.parameters)>numel(theta) % if the model parameter vector is longer than the loaded theta 
    nParameterDiff = numel(model.parameters)-numel(theta); % calulate the difference in parameter vector length 
    theta(end+1:end+nParameterDiff) = zeros(1,nParameterDiff); % add zeros to the unkown parameter values 
end

%% configure the total parameter estimation vector to contian different parameter sets for different experiments 
for i = 1:size(parameterMappingTable,1) % loop throught the parameter Mapping table.
    if any(parameterMappingTable.TotalParameterIndex{i}>numel(theta)) % if any of the TotalParameterIndex on th i:th row exceed the number of parameters in the loaded theta-vector. 
        theta(parameterMappingTable.TotalParameterIndex{i}) = theta(parameterMappingTable.ModelParameterIndex{i}); % set the inital value of any not previously existing parameter values to that of the corresponding parameter index .
        thetaOptParameterNames(parameterMappingTable.TotalParameterIndex{i}) = thetaOptParameterNames(parameterMappingTable.ModelParameterIndex{i}); % set the corresponding parameters names for the added parameters.
    end

    if any(parameterMappingTable.TotalParameterIndex{i}>numel(thetaOptParameterNames)) % if any of the TotalParameterIndex on th i:th row exceed the number of parameters to optimize 
        thetaOptParameterNames(parameterMappingTable.TotalParameterIndex{i}) = thetaOptParameterNames(parameterMappingTable.ModelParameterIndex{i}); % set the corresponding parameters names for the added parameters.
    end
end

%% Allow the user to specify the index of specific parameters to estimate. All other parameters will be fixated at the inital value.  
OptParamBool = false(1, numel(thetaOptParameterNames)); % preallocate a boolean value for every parameter. The boolean value is used to determine if a parameter is included in the parameter estimation or if it is keep constant. 1 = included, 0 = keept constant. 

maxNParametersInSetup = max(max(cell2mat(parameterMappingTable.TotalParameterIndex.')),numel(model.parameters)); % compare and get the larger numebr of the largest index in the totalParameter Idx or the number of parameters in the model.
OptParamIndexes = setdiff([1:maxNParametersInSetup],Index); % set index of the parameters to be included in parameter estimation. Default is the full maximal number of parameters in the current setup , depending on which experiments are choosen. 

OptParamBool(OptParamIndexes) = true; % set parameters to optimise to have a true flag

thetaOpt = theta(OptParamBool); % build the parameter vector that is passed to parameter estimation algorithm.
thetaConstants = theta(~OptParamBool); % build the parameter vector that is keept constant during parmater estimation.

thetaConstants(find(~OptParamBool)==Index) = referenceValue;
problem.x_0 = thetaOpt;

%% Generate parameter bounds
thetaOptParameterNames = thetaOptParameterNames(OptParamBool); % Names of parameters that are included in parameter estiamtion.
[lb,ub] = setLBandUB(thetaOptParameterNames);

problem.x_L       = lb; % essOPT uses a problem structure where crucial information is specified
problem.x_U       = ub;

%% MEIGO OPTIONS I (COMMON TO ALL SOLVERS):
opts.ndiverse   = 100;      %100; %500; %5; %
PErunTime = 1;

opts.maxtime    = PErunTime*3600;       % MAX-Time of optmization, i.e how long the optimization will last

opts.maxeval    = 1e8;      % max number of evals, i.e cost function calls
opts.log_var    = [];       % skip this

opts.local.solver = 'dhc';  % dhc'; %'fmincon'; %'nl2sol'; %'mix'; 
opts.local.finish = opts.local.solver; % uses the local solver to check the best p-vector
opts.local.bestx = 0;       
opts.local.balance = 0.5;   % how far from startguess the local search will push the params, 0-close, 1-far (Default 0.5)
opts.local.n1   = 1;        % Number of iterations before applying local search for the 1st time (Default 1)
opts.local.n2   = 2;        % Minimum number of iterations in the global phase between 2 local calls (Default 10) 

problem.f       = 'objectiveFunction'; % calls function that sets up the cost function call

% MEIGO OPTIONS II (FOR ESS AND MULTISTART):
opts.local.iterprint = 1; % prints what going on during optimization

% MEIGO OPTIONS III (FOR ESS ONLY):
opts.dim_refset   = 10; 

% OPTIONS AUTOMATICALLY SET AS A RESULT OF PREVIOUS OPTIONS:
opts.local.check_gradient_for_finish = 0; % DW: gradient checker
optim_algorithm = 'ess'; % 'multistart'; %  'cess'; 

if any(problem.x_0>ub) || any(problem.x_0<lb) % print a warning if any of the starting parameters are outside lb and ub
    warning(['Start value for parameters ',repmat('%d ',1,nnz(problem.x_0<lb)) ,'are outside the lower bound condition\n'],find(problem.x_0<lb))
    warning(['Start value for parameters ',repmat('%d ',1,nnz(problem.x_0>ub)) ,'are outside the upper bound condition\n'],find(problem.x_0>ub))
end

%% get the Chi-threshold
experimentFieldNames = fieldnames(DataStructure);
nDataPoints = 0;
for i = 1:numel(ExperimentIdx)
     for j=1:numel(variableNames{i})
         % For variable j in Experiment i add the number of data points for a time point > 0
        nDataPoints =  nDataPoints + nnz(DataStructure.(experimentFieldNames{ExperimentIdx(i)}).(variableNames{i}{j}).Time>0);
     end
end

warning('off','all') 
fprintf('#------------------------------------#\n\nRunning Parameter Estimation\nModel:\t%s\nData:\t%s\nRuntime set to:\t%d seconds\nInitiated on:\t%s\nParameter index:\t%d\nPolarity:\t%d\nStep:\t%d\nTheta reference Value:\t%0.4f\n#------------------------------------#\n',modelName, strjoin(saveString,' - '),opts.maxtime,datestr(now,'yymmdd_HHMMSS'),Index, polarity,step,referenceValue)
Results = MEIGO(problem,opts,optim_algorithm, modelName, thetaConstants, OptParamBool, circuitConstants, DataStructure, ExperimentIdx, variableNames ,modelObservablesIndex, parameterMappingTable, model,options, 0,[]);
warning('on','all') 

%% reconstruct full parameter vector
xbest = zeros(1,length(OptParamBool));
xbest(OptParamBool) = Results.xbest;
xbest(~OptParamBool) = thetaConstants;

theta = xbest;
%% get fval and simulationStructure of best solution
[fval,simulation] = objectiveFunction(theta(OptParamBool), modelName, thetaConstants, OptParamBool, circuitConstants, DataStructure, ExperimentIdx, variableNames,modelObservablesIndex, parameterMappingTable, model, options,0,[]);

%% save in designated location 
saveFileName = sprintf('Index_%d_polarity_%d_step_%d_fval_(%.2f)_%s_%s.mat',Index, polarity, step, fval,datestr(now,'yymmdd_HHMMSS'),modelName);

save([ResultLocation,saveFileName],'theta','fval','simulation','Results','parameterMappingTable','saveString','Index','referenceValue')
end

