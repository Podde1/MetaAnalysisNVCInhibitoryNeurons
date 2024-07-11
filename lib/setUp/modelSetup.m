function [circuitConstants, DataStructure, model, variableNames, saveString, modelObservablesIndex, parameterMappingTable ,options] = modelSetup(modelName,ExperimentIdx,variables2Fit)

circuitConstants = getCircuitConstants();

%% load data
load('CombinedDataStructure.mat', 'CombinedDataStructure')
DataStructure = CombinedDataStructure;

%% set up model 
model = getModelNames([modelName,'_syms']); % get the names of model states, parameters, constants, and observables from the syms file.

model.observables = strrep(model.observables,'_y',''); % remove the "_y" sufix from the model observable strings.

%% Set up the Experiment
% Experiment idecies for combined Data sturcute
%  1   {'Vazques2018_Excitatory_FL'         }
%  2   {'Vazques2018_Excitatory_2ms'        }     
%  3   {'Vazques2018_Excitatory_10ms'       }
%  4   {'Vazques2018_Excitatory_30ms'       }
%  5   {'Vazques2018_Inhibitory_FL'         }
%  6   {'Vazques2018_Inhibitory_2ms'        } 
%  7   {'Vazques2018_Inhibitory_10ms'       }
%  8   {'Vazques2018_Inhibitory_30ms'       }
%  9   {'Moon2021_Excitatory_1Hz_20s'       }
%  10  {'Moon2021_Excitatory_20Hz_20s'      }
%  11  {'Moon2021_Forepaw_4Hz_20s'          }
%  12  {'Moon2021_Inhibitory_1Hz_20s'       }
%  13  {'Moon2021_Inhibitory_20Hz_20s'      }
%  14  {'Moon2021_Inhibitory_1Hz_5s'        }
%  15  {'Moon2021_Inhibitory_20Hz_5s'       }
%  16  {'Lee2020_Whiskers_5Hz_2s'           }
%  17  {'Lee2020_Whiskers_5Hz_16s'          }
%  18  {'Lee2020_SST_20Hz_2s'               }
%  19  {'Lee2020_SST_20Hz_16s'              }
%  20  {'Lee2020_nNOS_20Hz_2s'              }
%  21  {'Lee2020_nNOS_20Hz_16s'             }

nExperiments = numel(ExperimentIdx); % get the number of experiments used

% prellocate variables 
StimulationParadigms = zeros(nExperiments,3); % the stimulation paradigms col1: frequency, col2: PulseWidth
variableNames = cell(1, nExperiments);
modelObservablesIndex = cell(1, nExperiments);
saveString = cell(1, nExperiments);

% Get the field names of the DataStructure
experimentFieldNames = fieldnames(DataStructure);

% Loop over each experiment index
for i = 1:nExperiments
    % Get the field names of the experiment-specific data structure
    variableFieldNames = fieldnames(DataStructure.(experimentFieldNames{ExperimentIdx(i)}));
    
    % Filter the variable field names based on variables2Fit
    variableNames{i} = variableFieldNames(contains(variableFieldNames, variables2Fit));
    
    % Loop over each variable name in the current experiment
    for j = 1:numel(variableNames{i})
        % Find the index of the variable name in the model observables
        modelObservablesIndex{i}(j) = find(strcmp(variableNames{i}(j), model.observables),1);
    end

    % Create the save string for the current experiments. Indicating which experiments are part of fit in the save file.
    saveString{i} = sprintf('%d_%s', ExperimentIdx(i), DataStructure.(experimentFieldNames{ExperimentIdx(i)}).meta.abbreviation);

    % Set a flag to differentiate between studies with different anesthesia. Moon and Lee has the same stimuli but different experimental settings and should not get the same parameters
    studyFlag = 0;
    if contains(experimentFieldNames{ExperimentIdx(i)}, 'Vazquez')
        studyFlag=1;
    elseif contains(experimentFieldNames{ExperimentIdx(i)}, 'Moon')
        studyFlag=2;
    elseif contains(experimentFieldNames{ExperimentIdx(i)}, 'Lee')
        studyFlag=3;
    end

    % Extract the stimulation parameters for the current experiment
    StimulationParadigms(i, :) = [DataStructure.(experimentFieldNames{ExperimentIdx(i)}).Stimulation.Frequency, DataStructure.(experimentFieldNames{ExperimentIdx(i)}).Stimulation.pulseWidth, studyFlag];

end
saveString = [unique(vertcat(variableNames{:}),'stable');saveString.']; %Add on eadditional string to the front of the save string. This string contains the unique variables that are part of the fit. This is based on variableNames rather than variables2Fit since variables2Fit might contain variables that are in no experiments.

%% set up unique stimuation parameters for differetn stimulation conditions

parametersAffectedByStimIndex = find(~cellfun(@isempty,regexp(model.parameters,'k_u[0-9]'))); %get the index of the parameters affected by stimulation paramdigm i.e. parametes maching k_u[0-9]
parametersAffectedByAnesthesiaIndex = find(~cellfun(@isempty,regexp(model.parameters,'k_u[0-9]')),1,'last')+1:find(strcmp(model.parameters,'k_Ca'))-1;% get the index of the parameters affected by the anesthesia paramdigm i.e. the parameters listed after k_ui but before k_Ca

parametersAffectedByStim = model.parameters(parametersAffectedByStimIndex); % get the names of the parameters affected by the stimulation paradigm
parametersAffectedByAnesthesia = model.parameters(parametersAffectedByAnesthesiaIndex);% get the names of the parameters affected by hte anesthesia paradigm 

[~,~,stimulationParadigmsIndex] = unique(StimulationParadigms,'rows','stable'); % Each unique rows in the stimulation Pradigm matrix indicates an additional set of stimulation parameters in the total parameter vector
[~,~,AnesthesiaParadigmIndex] = unique(extractBefore(experimentFieldNames(ExperimentIdx),'_'),'stable'); % get the anesthesia index of the experiments in ExperimentIdx.

nModelParameters = numel(model.parameters); % get the number of parameters in the model structure 
parameterIndexMapping = cell(nExperiments,2); % prallocate the cellstructure that holds the mapping between the model parameter vector and the total parameter optimization vector 
parameterMappingTable = table(parameterIndexMapping(:,1),parameterIndexMapping(:,2),'VariableNames',{'ModelParameterIndex','TotalParameterIndex'},'RowNames',experimentFieldNames(ExperimentIdx));

% set the first row/experiment of the table to be the default parameter vector.
parameterMappingTable.ModelParameterIndex{1,:} = [parametersAffectedByStimIndex,parametersAffectedByAnesthesiaIndex];
parameterMappingTable.TotalParameterIndex{1,:} = [parametersAffectedByStimIndex,parametersAffectedByAnesthesiaIndex];

parametersAffectedByAnesthesiaTableIndex = [1:numel(parametersAffectedByAnesthesiaIndex)]+numel(parametersAffectedByStimIndex);

for i = 2:nExperiments %loop over the remaining experiments 
    parameterMappingTable.ModelParameterIndex{i} = [parametersAffectedByStimIndex,parametersAffectedByAnesthesiaIndex]; % for row i set the modelParameter Idex. This will be the same for all experiments/rows
    if ismember(stimulationParadigmsIndex(i),stimulationParadigmsIndex(1:i-1)) % if any of the previous experiments have the same stimulation paradigm as experiment i, used those indexes for the parameteres affected by the stimulation paradigm 
        parameterMappingTable.TotalParameterIndex{i}(parametersAffectedByStimIndex) = parameterMappingTable.TotalParameterIndex{find(ismember(stimulationParadigmsIndex(1:i-1),stimulationParadigmsIndex(i)),1,'first')}(parametersAffectedByStimIndex);
    else % if the i:th experiment represents a new stimulation paradigm. in the corresponding indecies add the correct number of variable indices, based on the largest value in TotalParameterIndex or the number of parameters in the model. 
        parameterMappingTable.TotalParameterIndex{i}(parametersAffectedByStimIndex) = (1:numel(parametersAffectedByStim))+ max([cell2mat(parameterMappingTable.TotalParameterIndex.'),nModelParameters]);
    end
    
    if ismember(AnesthesiaParadigmIndex(i),AnesthesiaParadigmIndex(1:i-1))% if any of the previous experiments have the same anesthesia paradigm as experiment i, used those indexes for the parameteres affected by the anesthesia paradigm 
        parameterMappingTable.TotalParameterIndex{i}(parametersAffectedByAnesthesiaTableIndex) = parameterMappingTable.TotalParameterIndex{find(ismember(AnesthesiaParadigmIndex(1:i-1),AnesthesiaParadigmIndex(i)),1,'first')}(parametersAffectedByAnesthesiaTableIndex);
    else % if the i:th experiment represents a new anesthesia paradigm. In the corresponding indecies add the correct number of variable indices, based on the largest value in TotalParameterIndex or the number of parameters in the model. 
        parameterMappingTable.TotalParameterIndex{i}(parametersAffectedByAnesthesiaTableIndex) = (1:numel(parametersAffectedByAnesthesia))+ max([cell2mat(parameterMappingTable.TotalParameterIndex.'),nModelParameters]);
    end
end

%% Set ami options
options = amioption('sensi',0,...
    'maxsteps',1e3);

% alter simulation tolerances, DAE solver can not handle the default values
options.atol = 1e-5;
options.rtol = 1e-6;

end