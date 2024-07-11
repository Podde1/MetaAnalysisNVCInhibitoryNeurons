function [fullSimulation] = simulateSquareStenNVCmodel(theta,time,modelName, stimulationInfo,circuitConstants,inputTable,modelPropertyNames,options)

%% split time vector 
stimulationLength = stimulationInfo.Length;
stimulationPeriod = (1/stimulationInfo.Frequency); % Number of pulses per sec
stimulationPulseWidth = stimulationInfo.pulseWidth;   % stimulationInfo.Frequency*stimulationInfo.pulseWidth

periodTimePoints =[0:stimulationPeriod:(stimulationLength-stimulationPeriod)]+[0;stimulationPulseWidth];
stimulationTime = unique([periodTimePoints(:).',time(time<=stimulationLength),0]);

postStimulationTime  = time(time>stimulationLength)- stimulationLength;% create the time vector for the post stimulation period. 

%% set up model simulation handle 
simulateModel = str2func(['simulate_',modelName]);
simulateModel_steadyState = str2func(['simulate_',modelName,'_steadyState']);

% initiate result structure 
fullSimulation.status = 1;

modelIputLogicIndex = contains(inputTable.InputVariable,modelPropertyNames.constants); %get the logical index of the input variables in inputTable.InputVariable that exist in the model file 
IntraNeuronalSignallingStatesIdx = find(startsWith(modelPropertyNames.states,{'Glut','V1'}));
IntraNeuronalSignallingStatesIdx = [IntraNeuronalSignallingStatesIdx(1)+1:IntraNeuronalSignallingStatesIdx(end)-1];

%  get the state indecies of the 'NOvsm','PGE2vsm','NPYvsm' states
if ~isempty(regexp(modelName,'v2[0-9]', 'once'))
    NeuronalvasoactiveEffectsStatesIdx=[15 13 17];
else
    NeuronalvasoactiveEffectsStatesIdx=[19 15 21];
end
%% SS
NOvsm0 = 0; PGE2vsm0 = 0;  NPYvsm0 = 0; Ca_start = 10;
modelSteadyStateConstants = [inputTable.('Steady-state')(modelIputLogicIndex).', NOvsm0, PGE2vsm0, NPYvsm0, Ca_start, circuitConstants,zeros(1,7)];

steadyState_simulation = simulateModel_steadyState([1,5e6],theta,modelSteadyStateConstants,[],options);

% if steadystate simulation return a negative status or has not reach a good enough steady state - return an error to indicate that the simulation failed.
if or(or(steadyState_simulation.status<0, any(steadyState_simulation.diagnosis.xdot > 1e-8)), any(steadyState_simulation.x(end,IntraNeuronalSignallingStatesIdx)<0))
    fullSimulation.status = 0;
    return
end

% assaign values to constants in the stimulation simulation
options.x0 = steadyState_simulation.x(end,:).';

%% set up stimulation simulation  
TE = 3*10^-3;      
B0 = 15.2;

HbO_0 = steadyState_simulation.y(end,1);   HbR_0 = steadyState_simulation.y(end,2);
SaO2_0 = steadyState_simulation.y(end,3);  ScO2_0 = steadyState_simulation.y(end,4);  SvO2_0 = steadyState_simulation.y(end,5);

simTime = zeros(size(stimulationTime.'));
simX = zeros(numel(simTime),size(steadyState_simulation.x,2));
simY = zeros(numel(simTime),length(modelPropertyNames.observables));
simStatus = zeros(size(periodTimePoints));

% Set variables for model constants. Since these are constants they can be set befor the loop as one for stimulation on and one for stimulation off.
modelConstants =      [inputTable.('Stimulation')(modelIputLogicIndex).'     ,steadyState_simulation.x(end,NeuronalvasoactiveEffectsStatesIdx), Ca_start, circuitConstants, HbO_0, HbR_0, SaO2_0, ScO2_0, SvO2_0, TE, B0];
modelConstants_rest = [inputTable.('Post-stimulation')(modelIputLogicIndex).',steadyState_simulation.x(end,NeuronalvasoactiveEffectsStatesIdx), Ca_start, circuitConstants, HbO_0, HbR_0, SaO2_0, ScO2_0, SvO2_0, TE, B0];


%% simulate pulses
nPeriods = size(periodTimePoints,2)-1;
for i=1:nPeriods %loop over stimulation period, trigger a pulse at every period i.e. stimulationFreq times per sec 
    
    stimulationTimeIdx = and(stimulationTime>=periodTimePoints(1,i), stimulationTime<=periodTimePoints(2,i));
    tStim = stimulationTime(stimulationTimeIdx);
    stimulation_simulation = simulateModel(tStim-tStim(1),theta,modelConstants,[],options);

    options.x0 = stimulation_simulation.x(end,:)';

    restTimeIdx = and(stimulationTime>=periodTimePoints(2,i), stimulationTime<=periodTimePoints(1,i+1));
    tRest = stimulationTime(restTimeIdx);
    rest_simulation = simulateModel(tRest-tRest(1),theta,modelConstants_rest,[],options);
    options.x0 = rest_simulation.x(end,:)';

    simX(stimulationTimeIdx,:) = stimulation_simulation.x;
    simY(stimulationTimeIdx,:) = stimulation_simulation.y;
    simX(restTimeIdx,:) = rest_simulation.x;
    simY(restTimeIdx,:) = rest_simulation.y;

    if any(simStatus(:)<0) % if either the stimulation or rest simulations fail return a negative status return an error to indicate that the simulation failed. and return the value to reduce unnecessary simulations if a falied state is reached
        fullSimulation.status = 0;
        return
    end
end

%% simulate final pulse and post stimulation simulation 

stimulationTimeIdx = and(stimulationTime>=periodTimePoints(1,end), stimulationTime<=periodTimePoints(2,end));
tStim = stimulationTime(stimulationTimeIdx);
stimulation_simulation = simulateModel(tStim-tStim(1),theta,modelConstants,[],options);
options.x0 = stimulation_simulation.x(end,:)';

restTimeIdx = and(stimulationTime>=periodTimePoints(2,end), stimulationTime<=stimulationLength);
tRest = stimulationTime(restTimeIdx);
rest_simulation = simulateModel(tRest-tRest(1),theta,modelConstants_rest,[],options);
options.x0 = rest_simulation.x(end,:)';

postStimulation_simulation = simulateModel(postStimulationTime,theta,modelConstants_rest,[],options);

simX(stimulationTimeIdx,:) = stimulation_simulation.x;
simY(stimulationTimeIdx,:) = stimulation_simulation.y;
simX(restTimeIdx,:) = rest_simulation.x;
simY(restTimeIdx,:) = rest_simulation.y;

%% Combine simulations
StimulationTimeIndex = ismember(stimulationTime,time);

% remove the elements corrisponding to the 0 and stimulationLenght time points in the time vector 
fullSimulation.t = [stimulationTime(StimulationTimeIndex).'  ; postStimulation_simulation.t + stimulationTime(end)];
fullSimulation.x = [simX(StimulationTimeIndex,:); postStimulation_simulation.x];
fullSimulation.y = [simY(StimulationTimeIndex,:); postStimulation_simulation.y];
fullSimulation.steadyState = steadyState_simulation.x(end,:);

% if any of the simulations return a negative status return an error to indicate that the simulation failed.
if any([stimulation_simulation.status, postStimulation_simulation.status]<0)
    fullSimulation.status = 0;
end
end