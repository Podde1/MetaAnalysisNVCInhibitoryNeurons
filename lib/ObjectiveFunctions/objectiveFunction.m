function [fval,simulation] =  objectiveFunction(theta, modelName, thetaConstants, OptParamBool, circuitConstants, Data,ExperimentIdx,variableNames,modelObservablesIndex, parameterMappingTable, modelPropertyNames, options,plotLogic,fid)
if size(theta,2)>1;theta = theta.';end % ensure theta is a column vector

if any(~OptParamBool)
    if size(thetaConstants,2)>1;thetaConstants = thetaConstants.';end % ensure thetaConstants is a column vector
    
    % reconstruct full theta vector
    temp = zeros(1,length(OptParamBool)).'; %create a temporary vector to reconstruc the full parameter vector 
    temp(OptParamBool) = theta; % add the parameters that are evaluated by the algorithm to the temp vector 
    temp(~OptParamBool) = thetaConstants;  % add the parameter values that ar ekeept constant. 
    theta = temp; 
end

%% Set up time vectors 
stimulationStartTime = 0; % set the stimulation start time. Used to ignore any data points prior to stimulation
squareResiduals = cell(1,sum(cellfun('size',variableNames,1))); %find the number of varaibles to fit to in order to preallocate the number of cells 

experimentNames = fieldnames(Data); % get the field names from the data structure
TotalCount = 1;
ExperimentCount = 1;
fval = 0;

for i = ExperimentIdx
    stimulationLength = Data.(experimentNames{i}).Stimulation.Length; %Used to split the simulation into stimuli and post stimuli simulaiton
    stimulationInfo = Data.(experimentNames{i}).Stimulation; %Used to split the simulation into stimuli and post stimuli simulaiton
    %% set up the time vector for simulating experimetn i
    simulationTime=[];
    for j = 1:numel(variableNames{ExperimentCount})  % loop over the variables to fit in Experiment i. 
        simulationTime = [simulationTime,Data.(experimentNames{i}).(variableNames{ExperimentCount}{j}).Time]; % gather all the different time vectors into one   
    end
    highResStimTime = 0;
    if highResStimTime
        simulationTime = unique([simulationTime(simulationTime>stimulationStartTime),0, [simulationTime(find(simulationTime>0,1,'first')):0.01:stimulationInfo.Length]]);
    else
        simulationTime = unique([simulationTime(simulationTime>stimulationStartTime),0]); % sort out duplicate values and sort the vector. Also sort out values <0. add 0 such that this time poind can be includend in plots
    end
    %% set  up the parameter vector for simulating experiment i
    % calculate the number of parameters in the model based on the entries of the parameterMapping table. i.e. the number of model parameters is excatly 1 lower thatn the smallset value in the TotalParameterIndex-column that is also not in the ModelParameterIndex-column. 
   nModelParameters = min(setdiff(cell2mat(parameterMappingTable.TotalParameterIndex.'),cell2mat(parameterMappingTable.ModelParameterIndex.')))-1;
   if isempty(nModelParameters); nModelParameters = numel(theta);end % if there is no set difference between ModelParameterIndex and the TotalparameterIndex, get the length of the model parameter vector from theta. 
   modelTheta = theta(1:nModelParameters);

    % assign the parameters to use for the model simulation based on the row of the parameterMappingTable that maches the i:tn experiment name.
    modelTheta(parameterMappingTable(experimentNames{i},:).ModelParameterIndex{1}) = theta(parameterMappingTable(experimentNames{i},:).TotalParameterIndex{1});

    %% run a single simulation for Experiment i 
    simulation.(experimentNames{i}) = simulateSquareStenNVCmodel(modelTheta,simulationTime, modelName, stimulationInfo, circuitConstants, Data.(experimentNames{i}).input, modelPropertyNames,options);

    if ~simulation.(experimentNames{i}).status %if the simulation fails. Return a very high cost. 
        fval = fval + 1e29;
        return % currently this objFunction returns as soon as any simulation fails. One option is to complete all simulations, save the status flag and then modify fval if any has failed. 
    else
        %% calculate the wheighted squared residuals
        for j = 1:numel(variableNames{ExperimentCount}) % loop over the differnt variables to find the weighted squared residuals
            variableNameStr = variableNames{ExperimentCount}{j}; % get the jth variable name as a string 
            dataPointsIdx = Data.(experimentNames{i}).(variableNameStr).Time>=stimulationStartTime; % get the logical index of the data points that are after the stimulaiton starts
            ExperimentJ_timeidx = ismember(simulationTime,Data.(experimentNames{i}).(variableNameStr).Time);% get the logical idx of which time points in the large time vector corresponds to the time points of variable j. 
    
            % calculate the SEM wighted squared residuals for all data points in variable j.
            squareResiduals{TotalCount} = ((Data.(experimentNames{i}).(variableNameStr).Mean(dataPointsIdx).' - simulation.(experimentNames{i}).y(ExperimentJ_timeidx,modelObservablesIndex{ExperimentCount}(j)))./Data.(experimentNames{i}).(variableNameStr).SEM(dataPointsIdx).').^2;
      
          %% if selcted. plot the current fit
            if plotLogic %if the plot logic applies
                figure() % open a  new figure
                hold on
                % plot the datapoints of variable j in Experiment i as error bars 
                errorbar(Data.(experimentNames{i}).(variableNameStr).Time,Data.(experimentNames{i}).(variableNameStr).Mean,Data.(experimentNames{i}).(variableNameStr).SEM,' r*','LineWidth',2)
        
                plotTimeidx = or(ExperimentJ_timeidx.',simulation.(experimentNames{i}).t==0); % ensure that t=0 is included in simulation plot
                plot(simulation.(experimentNames{i}).t(plotTimeidx),simulation.(experimentNames{i}).y(plotTimeidx,modelObservablesIndex{ExperimentCount}(j)),'b-','LineWidth',2) % Plot the simulaiton of variable j in experiment i. 
                ax = gca;% get the current axes
                line([0, stimulationLength],repmat(ax.YLim(1),1,2),'LineWidth',5,'Color',[0 0 0])% plot a line at the bottom of the axes to represent the stimulation time
                title(strrep(strjoin([experimentNames(i);variableNames{ExperimentCount}(j)],' '),'_',' '))  % set the axes titel to the joint string of the experietn name and the variable name 
                hold off
            end
            
            TotalCount = TotalCount + 1;
        end
    end
    ExperimentCount = ExperimentCount + 1;
end

fval = sum(cell2mat(squareResiduals.')); % sum all squared residuals and the adhoc penalty. 

%% save acceptable parameters (used in MCMC sampling)
nDataPoints = numel(cell2mat(squareResiduals.')); % calculate the number of data points used for estimation
acceptCondition = fval <= chi2inv(0.95,nDataPoints);

if ((nargin >= 13) && ~isempty(fid)) && acceptCondition % set threshold at the total number of data points to gather as much as possible. stricter thresholds can be added at a later state
   fprintf(fid,'%4.10f %10.10f ',[fval, theta.']); fprintf(fid,'\n');
end

end
